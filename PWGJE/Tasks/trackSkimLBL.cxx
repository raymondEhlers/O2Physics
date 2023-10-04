// Copyright 2019-2023 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

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

// Track skimming task for Berkeley folks
//
// Author: Raymond Ehlers <raymond.ehlers@cern.ch>, LBL/UC Berkeley

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// Define track skim table
namespace o2::aod
{
namespace LBLSkim
{
// NOTE: We don't want a Collision Index ID, as this will require the full input AO2D, which we don't want!
DECLARE_SOA_INDEX_COLUMN(Collision, collision);
DECLARE_SOA_COLUMN(CollisionIndex, collisionIndex, int);
// Event level
// TODO: RJE: Check precision and ranges to ensure they're suitable for all values.
//            They could be tuned up or down based on requirements.
DECLARE_SOA_COLUMN(RunNumber, runNumber, int32_t);
DECLARE_SOA_COLUMN(Centrality, centrality, float32_t);
DECLARE_SOA_COLUMN(EventPlane, eventPlane, float32_t);
// Track / cluster level
DECLARE_SOA_COLUMN(Pt, pt, float32_t);
DECLARE_SOA_COLUMN(Eta, eta, float32_t);
DECLARE_SOA_COLUMN(Phi, phi, float32_t);
DECLARE_SOA_COLUMN(Energy, energy, float32_t);
DECLARE_SOA_COLUMN(ParticleID, particleID, int32_t);
DECLARE_SOA_COLUMN(Label, label, int32_t);
DECLARE_SOA_COLUMN(EncodedInformation, encodedInformation, int32_t);
} // namespace LBLSkim

// TODO: RJE: MC vs Data. We want to only include MC fields if they're meaningful. Probably best to define separate tables
DECLARE_SOA_TABLE(EventSkimLBL, "AOD", "EVTSKIMLBLV1",
                  o2::soa::Index<>,
                  LBLSkim::CollisionIndex,
                  LBLSkim::RunNumber
                  //LBLSkim::Centrality,
                  //LBLSkim::EventPlane
);
DECLARE_SOA_TABLE(TrackSkimLBL, "AOD", "TRKSKIMLBLV1",
                  o2::soa::Index<>,
                  LBLSkim::CollisionIndex,
                  LBLSkim::Pt,
                  LBLSkim::Eta,
                  LBLSkim::Phi,
                  // These fields should be conditional
                  LBLSkim::ParticleID,
                  LBLSkim::Label,
                  LBLSkim::EncodedInformation
                  );
                  // We can add whatever is desired here, including conditionally
DECLARE_SOA_TABLE(ClusterSkimLBL, "AOD", "CLUSSKIMLBLV1",
                  o2::soa::Index<>,
                  LBLSkim::CollisionIndex,
                  LBLSkim::Eta,
                  LBLSkim::Phi,
                  LBLSkim::Energy,
                  // These fields should be conditional
                  LBLSkim::Label,
                  LBLSkim::EncodedInformation
                  );
// For eg photon candidates
DECLARE_SOA_TABLE(PhotonSkimLBL, "AOD", "PHOSKIMLBLV1",
                  o2::soa::Index<>,
                  LBLSkim::CollisionIndex,
                  LBLSkim::Eta,
                  LBLSkim::Phi,
                  LBLSkim::Energy);
} // namespace o2::aod

using Tracks = o2::soa::Filtered<o2::soa::Join<o2::aod::Tracks, o2::aod::TrackSelection>>;
using Clusters = o2::soa::Filtered<o2::aod::EMCALClusters>;

struct TrackSkimLBL {
  Produces<o2::aod::TrackSkimLBL> trackSkim;
  Produces<o2::aod::ClusterSkimLBL> clusterSkim;
  Produces<o2::aod::PhotonSkimLBL> photonSkim;

  // Skim level configurables
  Configurable<bool> includeTracks{"includeTracks", true, "Whether to include tracks in the skim"};
  Configurable<bool> includeClusters{"includeClusters", false, "Whether to include cluster in the skim"};

  // Event level configurables
  Configurable<std::string> eventSelections{"eventSelections", "sel8", "choose event selection"};
  Configurable<float> vertexZCut{"vertexZCut", 10.0f, "Accepted z-vertex range"};

  // Track level configurables
  Configurable<std::string> trackSelections{"trackSelections", "globalTracks", "set track selections"};
  Configurable<std::string> particleSelections{"particleSelections", "PhysicalPrimary", "set particle selections"};
  Configurable<float> trackPtMin{"trackPtMin", 0.15, "minimum track pT"};
  Configurable<float> trackPtMax{"trackPtMax", 1000.0, "maximum track pT"};
  Configurable<float> trackEtaMin{"trackEtaMin", -0.9, "minimum track eta"};
  Configurable<float> trackEtaMax{"trackEtaMax", 0.9, "maximum track eta"};
  Configurable<float> trackPhiMin{"trackPhiMin", -999, "minimum track phi"};
  Configurable<float> trackPhiMax{"trackPhiMax", 999, "maximum track phi"};

  // Cluster level observables
  Configurable<std::string> clusterDefinitionS{"clusterDefinition", "kV3Default", "cluster definition to be selected, e.g. V3Default"};
  Configurable<float> clusterEtaMin{"clusterEtaMin", -0.7, "minimum cluster eta"}; // For ECMAL: |eta| < 0.7, phi = 1.40 - 3.26
  Configurable<float> clusterEtaMax{"clusterEtaMax", 0.7, "maximum cluster eta"};  // For ECMAL: |eta| < 0.7, phi = 1.40 - 3.26
  Configurable<float> clusterPhiMin{"clusterPhiMin", -999, "minimum cluster phi"};
  Configurable<float> clusterPhiMax{"clusterPhiMax", 999, "maximum cluster phi"};
  Configurable<float> clusterEnergyMin{"clusterEnergyMin", 0.5, "minimum cluster energy in EMCAL (GeV)"};
  Configurable<float> clusterTimeMin{"clusterTimeMin", -999., "minimum Cluster time (ns)"};
  Configurable<float> clusterTimeMax{"clusterTimeMax", 999., "maximum Cluster time (ns)"};
  Configurable<bool> clusterRejectExotics{"clusterRejectExotics", true, "Reject exotic clusters"};

  std::string trackSelection;
  std::string eventSelection;
  std::string particleSelection;

  void init(InitContext const&)
  {
    // TODO: Still required to cast configurables into members?
    trackSelection = static_cast<std::string>(trackSelections);
    eventSelection = static_cast<std::string>(eventSelections);
    particleSelection = static_cast<std::string>(particleSelections);
  }

  // Filters
  o2::aod::EMCALClusterDefinition clusterDefinition = o2::aod::emcalcluster::getClusterDefinitionFromString(clusterDefinitionS.value);
  Filter collisionFilter = nabs(aod::collision::posZ) < vertexZCut;
  Filter trackCuts = (aod::track::pt >= trackPtMin && aod::track::pt < trackPtMax && aod::track::eta > trackEtaMin && aod::track::eta < trackEtaMax && aod::track::phi >= trackPhiMin && aod::track::phi <= trackPhiMax); // do we need eta cut both here and in the global selection?
  Filter partCuts = (aod::mcparticle::pt >= trackPtMin && aod::mcparticle::pt < trackPtMax && aod::mcparticle::eta > trackEtaMin && aod::mcparticle::eta < trackEtaMax);
  Filter clusterFilter = (o2::aod::emcalcluster::definition == static_cast<int>(clusterDefinition) && aod::emcalcluster::eta > clusterEtaMin && aod::emcalcluster::eta < clusterEtaMax && aod::emcalcluster::phi >= clusterPhiMin && aod::emcalcluster::phi <= clusterPhiMax && aod::emcalcluster::energy >= clusterEnergyMin && aod::emcalcluster::time > clusterTimeMin && aod::emcalcluster::time < clusterTimeMax && (clusterRejectExotics && aod::emcalcluster::isExotic != true));

  // Process functions
  void processTracks(
    soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>::iterator const& collision,
    Tracks const& tracks
  )
  {
    if (!selectCollision(collision, eventSelection)) {
      return;
    }

    for (auto& track : tracks) {
      trackSkim(collision.globalIndex(), track.pt(), track.eta(), track.phi());
    }
  }
  PROCESS_SWITCH(TrackSkimLBL, processTracks, "Tracks", true);

  void processClusters(
    soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>::iterator const& collision,
    JetClusters const& clusters
  )
  {
    if (!selectCollision(collision, eventSelection)) {
      return;
    }

    for (auto& cluster : clusters) {
      // TODO: RJE: Calculate properly with vertex correction
      float energy = std::sqrt(cluster.p() * cluster.p() + mPionSquared);
      clusterSkim(collision.globalIndex(), cluster.eta(), cluster.phi(), energy);
    }
  }
  PROCESS_SWITCH(TrackSkimLBL, processClusters, "Clusters", true);

  // TODO: RJE: Photon (pre)selection? Or just build from clusters?
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<TrackSkimLBL>(cfgc, TaskName{"track-skim-LBL"})};
}