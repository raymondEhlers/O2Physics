// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file candidateSelectorDsToKKPi.cxx
/// \brief Ds± → K± K∓ π± selection task
///
/// \author Fabio Catalano <fabio.catalano@cern.ch>, Universita and INFN Torino
/// \author Stefano Politano <stefano.politano@cern.ch>, Politecnico and INFN Torino

#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include "Common/Core/TrackSelectorPID.h"

#include "PWGHF/Core/HfMlResponse.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::aod::hf_cand_3prong;
using namespace o2::analysis::hf_cuts_ds_to_k_k_pi;

/// Struct for applying Ds to KKpi selection cuts
struct HfCandidateSelectorDsToKKPi {
  Produces<aod::HfSelDsToKKPi> hfSelDsToKKPiCandidate;
  Produces<aod::HfMlDsToKKPi> hfMlDsToKKPiCandidate;

  Configurable<double> ptCandMin{"ptCandMin", 1., "Lower bound of candidate pT"};
  Configurable<double> ptCandMax{"ptCandMax", 36., "Upper bound of candidate pT"};
  // TPC PID
  Configurable<double> ptPidTpcMin{"ptPidTpcMin", 0.15, "Lower bound of track pT for TPC PID"};
  Configurable<double> ptPidTpcMax{"ptPidTpcMax", 20., "Upper bound of track pT for TPC PID"};
  Configurable<double> nSigmaTpcMax{"nSigmaTpcMax", 3., "Nsigma cut on TPC"};
  //  TOF PID
  Configurable<double> ptPidTofMin{"ptPidTofMin", 0.15, "Lower bound of track pT for TOF PID"};
  Configurable<double> ptPidTofMax{"ptPidTofMax", 20., "Upper bound of track pT for TOF PID"};
  Configurable<double> nSigmaTofMax{"nSigmaTofMax", 3., "Nsigma cut on TOF"};
  // topological cuts
  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{hf_cuts_ds_to_k_k_pi::vecBinsPt}, "pT bin limits"};
  Configurable<LabeledArray<double>> cuts{"cuts", {hf_cuts_ds_to_k_k_pi::cuts[0], nBinsPt, nCutVars, labelsPt, labelsCutVar}, "Ds candidate selection per pT bin"};
  // QA switch
  Configurable<bool> activateQA{"activateQA", false, "Flag to enable QA histogram"};
  // ML inference
  Configurable<bool> applyMl{"applyMl", false, "Flag to apply ML selections"};
  Configurable<std::vector<double>> binsPtMl{"binsPtMl", std::vector<double>{hf_cuts_ml::vecBinsPt}, "pT bin limits for ML application"};
  Configurable<std::vector<int>> cutDirMl{"cutDirMl", std::vector<int>{hf_cuts_ml::vecCutDir}, "Whether to reject score values greater or smaller than the threshold"};
  Configurable<LabeledArray<double>> cutsMl{"cutsMl", {hf_cuts_ml::cuts[0], hf_cuts_ml::nBinsPt, hf_cuts_ml::nCutScores, hf_cuts_ml::labelsPt, hf_cuts_ml::labelsCutScore}, "ML selections per pT bin"};
  Configurable<int8_t> nClassesMl{"nClassesMl", (int8_t)hf_cuts_ml::nCutScores, "Number of classes in ML model"};
  // CCDB configuration
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> modelPathsCCDB{"modelPathsCCDB", "EventFiltering/PWGHF/BDTD0", "Path on CCDB"};
  Configurable<std::vector<std::string>> onnxFileNames{"onnxFileNames", std::vector<std::string>{"ModelHandler_onnx_D0ToKPi.onnx"}, "ONNX file names for each pT bin (if not from CCDB full path)"};
  Configurable<int64_t> timestampCCDB{"timestampCCDB", -1, "timestamp of the ONNX file for ML model used to query in CCDB"};
  Configurable<bool> loadModelsFromCCDB{"loadModelsFromCCDB", false, "Flag to enable or disable the loading of models from CCDB"};

  o2::analysis::HfMlResponse<float> hfMlResponse;
  std::vector<float> outputMl = {};

  o2::ccdb::CcdbApi ccdbApi;

  TrackSelectorPi selectorPion;
  TrackSelectorKa selectorKaon;

  using TracksSel = soa::Join<aod::TracksWDca, aod::TracksPidPi, aod::TracksPidKa>;

  HistogramRegistry registry{"registry"};

  void init(InitContext const&)
  {
    selectorPion.setRangePtTpc(ptPidTpcMin, ptPidTpcMax);
    selectorPion.setRangeNSigmaTpc(-nSigmaTpcMax, nSigmaTpcMax);
    selectorPion.setRangePtTof(ptPidTofMin, ptPidTofMax);
    selectorPion.setRangeNSigmaTof(-nSigmaTofMax, nSigmaTofMax);
    selectorKaon = selectorPion;

    if (activateQA) {
      constexpr int kNBinsSelections = 1 + aod::SelectionStep::NSelectionSteps;
      std::string labels[kNBinsSelections];
      labels[0] = "No selection";
      labels[1 + aod::SelectionStep::RecoSkims] = "Skims selection";
      labels[1 + aod::SelectionStep::RecoTopol] = "Skims & Topological selections";
      labels[1 + aod::SelectionStep::RecoPID] = "Skims & Topological & PID selections";
      labels[1 + aod::SelectionStep::RecoMl] = "ML selection";
      static const AxisSpec axisSelections = {kNBinsSelections, 0.5, kNBinsSelections + 0.5, ""};
      registry.add("hSelections", "Selections;;#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {axisSelections, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
      for (int iBin = 0; iBin < kNBinsSelections; ++iBin) {
        registry.get<TH2>(HIST("hSelections"))->GetXaxis()->SetBinLabel(iBin + 1, labels[iBin].data());
      }
    }

    if (applyMl) {
      hfMlResponse.configure(binsPtMl, cutsMl, cutDirMl, nClassesMl);
      if (loadModelsFromCCDB) {
        ccdbApi.init(ccdbUrl);
        hfMlResponse.setModelPathsCCDB(onnxFileNames, ccdbApi, modelPathsCCDB.value, timestampCCDB);
      } else {
        hfMlResponse.setModelPathsLocal(onnxFileNames);
      }
      hfMlResponse.init();
      outputMl.assign(((std::vector<int>)cutDirMl).size(), -1.f); // dummy value for ML output
    }
  }

  /// Candidate selections independent from the daugther-mass hypothesis
  /// \param candidate is candidate
  /// \return true if candidate passes all cuts
  template <typename T1>
  bool selection(const T1& candidate)
  {
    auto candpT = candidate.pt();
    int pTBin = findBin(binsPt, candpT);
    if (pTBin == -1) {
      return false;
    }

    if (candpT < ptCandMin || candpT > ptCandMax) { // check that the candidate pT is within the analysis range
      return false;
    }
    if (candidate.decayLength() < cuts->get(pTBin, "decay length")) {
      return false;
    }
    if (candidate.decayLengthXYNormalised() < cuts->get(pTBin, "normalized decay length XY")) {
      return false;
    }
    if (candidate.cpa() < cuts->get(pTBin, "cos pointing angle")) {
      return false;
    }
    if (candidate.cpaXY() < cuts->get(pTBin, "cos pointing angle XY")) {
      return false;
    }
    if (std::abs(candidate.impactParameterXY()) > cuts->get(pTBin, "impact parameter XY")) {
      return false;
    }
    return true;
  }

  /// Candidate selections for the KKPi daugther-mass hypothesis
  /// \param candidate is candidate
  /// \param trackKaon1 is the first track with the kaon hypothesis
  /// \param trackKaon2 is the second track with the kaon hypothesis
  /// \param trackPion is the track with the pion hypothesis
  /// \return true if candidate passes all cuts
  template <typename T1, typename T2>
  bool selectionKKPi(const T1& candidate, const T2& trackKaon1, const T2& trackKaon2, const T2& trackPion)
  {
    int pTBin = findBin(binsPt, candidate.pt());
    if (pTBin == -1) {
      return false;
    }

    if (trackKaon1.pt() < cuts->get(pTBin, "pT K") || trackKaon2.pt() < cuts->get(pTBin, "pT K") || trackPion.pt() < cuts->get(pTBin, "pT Pi")) {
      return false;
    }
    if (std::abs(invMassDsToKKPi(candidate) - RecoDecay::getMassPDG(pdg::Code::kDS)) > cuts->get(pTBin, "deltaM")) {
      return false;
    }
    if (deltaMassPhiDsToKKPi(candidate) > cuts->get(pTBin, "deltaM Phi")) {
      return false;
    }
    if (std::abs(cos3PiKDsToKKPi(candidate)) < cuts->get(pTBin, "cos^3 theta_PiK")) {
      return false;
    }
    return true;
  }

  /// Candidate selections for the PiKK daugther-mass hypothesis
  /// \param candidate is candidate
  /// \param trackPion is the track with the pion hypothesis
  /// \param trackKaon1 is the first track with the kaon hypothesis
  /// \param trackKaon2 is the second track with the kaon hypothesis
  /// \return true if candidate passes all cuts
  template <typename T1, typename T2>
  bool selectionPiKK(const T1& candidate, const T2& trackPion, const T2& trackKaon1, const T2& trackKaon2)
  {
    int pTBin = findBin(binsPt, candidate.pt());
    if (pTBin == -1) {
      return false;
    }

    if (trackKaon1.pt() < cuts->get(pTBin, "pT K") || trackKaon2.pt() < cuts->get(pTBin, "pT K") || trackPion.pt() < cuts->get(pTBin, "pT Pi")) {
      return false;
    }
    if (std::abs(invMassDsToPiKK(candidate) - RecoDecay::getMassPDG(pdg::Code::kDS)) > cuts->get(pTBin, "deltaM")) {
      return false;
    }
    if (deltaMassPhiDsToPiKK(candidate) > cuts->get(pTBin, "deltaM Phi")) {
      return false;
    }
    if (std::abs(cos3PiKDsToPiKK(candidate)) < cuts->get(pTBin, "cos^3 theta_PiK")) {
      return false;
    }
    return true;
  }

  void process(aod::HfCand3Prong const& candidates,
               TracksSel const&)
  {
    // looping over 3-prong candidates
    for (const auto& candidate : candidates) {

      // final selection flag:
      auto statusDsToKKPi = 0;
      auto statusDsToPiKK = 0;

      if (!(candidate.hfflag() & 1 << DecayType::DsToKKPi)) {
        hfSelDsToKKPiCandidate(statusDsToKKPi, statusDsToPiKK);
        if (applyMl) {
          hfMlDsToKKPiCandidate(outputMl);
        }
        if (activateQA) {
          registry.fill(HIST("hSelections"), 1, candidate.pt());
        }
        continue;
      }
      SETBIT(statusDsToKKPi, aod::SelectionStep::RecoSkims);
      SETBIT(statusDsToPiKK, aod::SelectionStep::RecoSkims);
      if (activateQA) {
        registry.fill(HIST("hSelections"), 2 + aod::SelectionStep::RecoSkims, candidate.pt());
      }

      auto trackPos1 = candidate.prong0_as<TracksSel>(); // positive daughter (negative for the antiparticles)
      auto trackNeg = candidate.prong1_as<TracksSel>();  // negative daughter (positive for the antiparticles)
      auto trackPos2 = candidate.prong2_as<TracksSel>(); // positive daughter (negative for the antiparticles)

      // topological selections
      if (!selection(candidate)) {
        hfSelDsToKKPiCandidate(statusDsToKKPi, statusDsToPiKK);
        if (applyMl) {
          hfMlDsToKKPiCandidate(outputMl);
        }
        continue;
      }

      bool topolDsToKKPi = selectionKKPi(candidate, trackPos1, trackNeg, trackPos2);
      bool topolDsToPiKK = selectionPiKK(candidate, trackPos1, trackNeg, trackPos2);
      if (!topolDsToKKPi && !topolDsToPiKK) {
        hfSelDsToKKPiCandidate(statusDsToKKPi, statusDsToPiKK);
        if (applyMl) {
          hfMlDsToKKPiCandidate(outputMl);
        }
        continue;
      }
      if (topolDsToKKPi) {
        SETBIT(statusDsToKKPi, aod::SelectionStep::RecoTopol);
      }
      if (topolDsToPiKK) {
        SETBIT(statusDsToPiKK, aod::SelectionStep::RecoTopol);
      }
      if (activateQA) {
        registry.fill(HIST("hSelections"), 2 + aod::SelectionStep::RecoTopol, candidate.pt());
      }

      // track-level PID selection
      int pidTrackPos1Pion = selectorPion.statusTpcOrTof(trackPos1);
      int pidTrackPos1Kaon = selectorKaon.statusTpcOrTof(trackPos1);
      int pidTrackPos2Pion = selectorPion.statusTpcOrTof(trackPos2);
      int pidTrackPos2Kaon = selectorKaon.statusTpcOrTof(trackPos2);
      int pidTrackNegKaon = selectorKaon.statusTpcOrTof(trackNeg);

      bool pidDsToKKPi = !(pidTrackPos1Kaon == TrackSelectorPID::Rejected ||
                           pidTrackNegKaon == TrackSelectorPID::Rejected ||
                           pidTrackPos2Pion == TrackSelectorPID::Rejected);

      bool pidDsToPiKK = !(pidTrackPos1Pion == TrackSelectorPID::Rejected ||
                           pidTrackNegKaon == TrackSelectorPID::Rejected ||
                           pidTrackPos2Kaon == TrackSelectorPID::Rejected);

      if (!pidDsToKKPi && !pidDsToPiKK) {
        hfSelDsToKKPiCandidate(statusDsToKKPi, statusDsToPiKK);
        if (applyMl) {
          hfMlDsToKKPiCandidate(outputMl);
        }
        continue;
      }
      if (pidDsToKKPi) {
        SETBIT(statusDsToKKPi, aod::SelectionStep::RecoPID);
      }
      if (pidDsToPiKK) {
        SETBIT(statusDsToPiKK, aod::SelectionStep::RecoPID);
      }
      if (activateQA) {
        registry.fill(HIST("hSelections"), 2 + aod::SelectionStep::RecoPID, candidate.pt());
      }

      if (applyMl) {
        // ML selections
        std::vector<float> inputFeatures{candidate.ptProng0(),
                                         trackPos1.dcaXY(),
                                         trackPos1.dcaZ(),
                                         candidate.ptProng1(),
                                         trackNeg.dcaXY(),
                                         trackNeg.dcaZ(),
                                         candidate.ptProng2(),
                                         trackPos2.dcaXY(),
                                         trackPos2.dcaZ()};

        bool isSelectedMl = hfMlResponse.isSelectedMl(inputFeatures, candidate.pt(), outputMl);
        hfMlDsToKKPiCandidate(outputMl);

        if (!isSelectedMl) {
          hfSelDsToKKPiCandidate(statusDsToKKPi, statusDsToPiKK);
          continue;
        }
        if (pidDsToKKPi) {
          SETBIT(statusDsToKKPi, aod::SelectionStep::RecoMl);
        }
        if (pidDsToPiKK) {
          SETBIT(statusDsToPiKK, aod::SelectionStep::RecoMl);
        }
        if (activateQA) {
          registry.fill(HIST("hSelections"), 2 + aod::SelectionStep::RecoMl, candidate.pt());
        }
      }

      hfSelDsToKKPiCandidate(statusDsToKKPi, statusDsToPiKK);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfCandidateSelectorDsToKKPi>(cfgc)};
}
