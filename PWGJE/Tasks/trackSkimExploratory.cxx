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

// Exploratory (slim) skim for testing viable approaches for a track skim
//
// Author: Raymond Ehlers <raymond.ehlers@cern.ch>, LBL/UC Berkeley

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"

#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/DataModel/JetReducedData.h"
//#include "PWGJE/TableProducer/jetfinder.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// TODO: Format!

// Define track skim table
namespace o2::aod
{
namespace ExploratorySkim
{
// NOTE: We don't want a Collision Index ID, as this will require the full input AO2D, which we don't want!
DECLARE_SOA_INDEX_COLUMN(Collision, collision);
DECLARE_SOA_COLUMN(CollisionIndex, collisionIndex, int64_t);
// Event level
// TODO: RJE: Check precision and ranges to ensure they're suitable for all values.
//            They could be tuned up or down based on requirements.
// TODO: RJE: Removed run number since it requires BC, which makes the rest more difficult. Can add back later
// DECLARE_SOA_COLUMN(RunNumber, runNumber, int32_t);
DECLARE_SOA_COLUMN(Centrality, centrality, float);
DECLARE_SOA_COLUMN(EventPlane, eventPlane, float);
DECLARE_SOA_COLUMN(EventWeight, eventWeight, float);
// Track / cluster level
namespace Rec
{
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(Phi, phi, float);
DECLARE_SOA_COLUMN(Energy, energy, float);
} // namespace Rec
namespace Part
{
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(Phi, phi, float);
DECLARE_SOA_COLUMN(Energy, energy, float);
} // namespace Part
DECLARE_SOA_COLUMN(ParticleID, particleID, int32_t);
DECLARE_SOA_COLUMN(Label, label, int32_t);
DECLARE_SOA_COLUMN(EncodedInformation, encodedInformation, int32_t);
} // namespace ExploratorySkim

/**********************************
 * Define the skim tables from here.
 *
 * NOTE: RJE: I manually defined tables for now, but for better future consistency, we might consider definining them
 *            with macros. That way, we can ensure the format will match up (since there's no concept of table inheritance)
 *            while still only including the MC fields conditionally.
 ***********************************/
/////////////////
// Event level tables
/////////////////
DECLARE_SOA_TABLE(EventSkimExploratory, "AOD", "EVTSKIMV1",
                  o2::soa::Index<>,
                  ExploratorySkim::CollisionIndex
                  // ExploratorySkim::RunNumber
                  // ExploratorySkim::Centrality,
                  // ExploratorySkim::EventPlane
);
DECLARE_SOA_TABLE(MCEventSkimExploratory, "AOD", "MCEVTSKIMV1",
                  o2::soa::Index<>,
                  ExploratorySkim::CollisionIndex,
                  // ExploratorySkim::RunNumber,
                  // ExploratorySkim::Centrality,
                  // ExploratorySkim::EventPlane
                  //  MC only fields
                  ExploratorySkim::EventWeight);
/////////////////
// Track tables
/////////////////
DECLARE_SOA_TABLE(TrackSkimExploratory, "AOD", "TRKSKIMV1",
                  o2::soa::Index<>,
                  ExploratorySkim::CollisionIndex,
                  ExploratorySkim::Rec::Pt,
                  ExploratorySkim::Rec::Eta,
                  ExploratorySkim::Rec::Phi);
// MC det and part level
// TODO: RJE: Need to be able to turn the encoded info field on and off
DECLARE_SOA_TABLE(MCDetLevelTrackSkimExploratory, "AOD", "MCDTRKSKIMV1",
                  o2::soa::Index<>,
                  ExploratorySkim::CollisionIndex,
                  ExploratorySkim::Rec::Pt,
                  ExploratorySkim::Rec::Eta,
                  ExploratorySkim::Rec::Phi,
                  // MC only fields
                  ExploratorySkim::ParticleID,
                  ExploratorySkim::Label,
                  ExploratorySkim::EncodedInformation);
// Same as det level, but with a separate name to ensure the tables don't conflict
DECLARE_SOA_TABLE(MCPartLevelTrackSkimExploratory, "AOD", "MCPTRKSKIMV1",
                  o2::soa::Index<>,
                  ExploratorySkim::CollisionIndex,
                  ExploratorySkim::Part::Pt,
                  ExploratorySkim::Part::Eta,
                  ExploratorySkim::Part::Phi,
                  // MC only fields
                  ExploratorySkim::ParticleID,
                  ExploratorySkim::Label,
                  ExploratorySkim::EncodedInformation);
/////////////////
// Cluster tables
/////////////////
DECLARE_SOA_TABLE(ClusterSkimExploratory, "AOD", "CLUSSKIMV1",
                  o2::soa::Index<>,
                  ExploratorySkim::CollisionIndex,
                  ExploratorySkim::Rec::Eta,
                  ExploratorySkim::Rec::Phi,
                  ExploratorySkim::Rec::Energy);
DECLARE_SOA_TABLE(MCDetLevelClusterSkimExploratory, "AOD", "MCDCLUSSKIMV1",
                  o2::soa::Index<>,
                  ExploratorySkim::CollisionIndex,
                  ExploratorySkim::Rec::Eta,
                  ExploratorySkim::Rec::Phi,
                  ExploratorySkim::Rec::Energy,
                  // MC only fields
                  ExploratorySkim::Label,
                  ExploratorySkim::EncodedInformation);
} // namespace o2::aod

using Tracks = o2::soa::Filtered<o2::soa::Join<o2::aod::Tracks, o2::aod::TrackSelection>>;
using Clusters = o2::soa::Filtered<o2::aod::EMCALClusters>;

struct TrackSkimExploratory {
  // We create all tables here, and then we'll conditionally fill them based on the input parameters
  Produces<o2::aod::EventSkimExploratory> eventSkim;
  Produces<o2::aod::TrackSkimExploratory> trackSkim;
  Produces<o2::aod::ClusterSkimExploratory> clusterSkim;
  Produces<o2::aod::MCEventSkimExploratory> eventSkimMC;
  Produces<o2::aod::MCDetLevelTrackSkimExploratory> trackSkimMCD;
  Produces<o2::aod::MCPartLevelTrackSkimExploratory> trackSkimMCP;
  Produces<o2::aod::MCDetLevelClusterSkimExploratory> clusterSkimMCD;

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
  Configurable<float> clusterEtaMin{"clusterEtaMin", -0.7, "minimum cluster eta"}; // For EMCAL: |eta| < 0.7, phi = 1.40 - 3.26
  Configurable<float> clusterEtaMax{"clusterEtaMax", 0.7, "maximum cluster eta"};  // For EMCAL: |eta| < 0.7, phi = 1.40 - 3.26
  Configurable<float> clusterPhiMin{"clusterPhiMin", -999, "minimum cluster phi"};
  Configurable<float> clusterPhiMax{"clusterPhiMax", 999, "maximum cluster phi"};
  Configurable<float> clusterEnergyMin{"clusterEnergyMin", 0.5, "minimum cluster energy in EMCAL (GeV)"};
  Configurable<float> clusterTimeMin{"clusterTimeMin", -999., "minimum Cluster time (ns)"};
  Configurable<float> clusterTimeMax{"clusterTimeMax", 999., "maximum Cluster time (ns)"};
  Configurable<bool> clusterRejectExotics{"clusterRejectExotics", true, "Reject exotic clusters"};

  // Track and event selection
  std::string eventSelection;
  std::string trackSelection;
  std::string particleSelection;
  int eventSelection = -1;
  int trackSelection = -1;
  int particleSelection = -1;

  void init(InitContext const&)
  {
    eventSelection = JetDerivedDataUtilities::initialiseEventSelection(static_cast<std::string>(eventSelections));
    trackSelection = JetDerivedDataUtilities::initialiseTrackSelection(static_cast<std::string>(trackSelections));
    particleSelection = JetDerivedDataUtilities::initialiseTrackSelection(static_cast<std::string>(particleSelections));
  }

  // Filters
  o2::aod::EMCALClusterDefinition clusterDefinition = o2::aod::emcalcluster::getClusterDefinitionFromString(clusterDefinitionS.value);
  Filter collisionFilter = nabs(aod::collision::posZ) < vertexZCut;
  Filter trackCuts = (aod::track::pt >= trackPtMin && aod::track::pt < trackPtMax && aod::track::eta > trackEtaMin && aod::track::eta < trackEtaMax && aod::track::phi >= trackPhiMin && aod::track::phi <= trackPhiMax); // do we need eta cut both here and in the global selection?
  Filter partCuts = (aod::mcparticle::pt >= trackPtMin && aod::mcparticle::pt < trackPtMax && aod::mcparticle::eta > trackEtaMin && aod::mcparticle::eta < trackEtaMax);
  Filter clusterFilter = (o2::aod::emcalcluster::definition == static_cast<int>(clusterDefinition) && aod::emcalcluster::eta > clusterEtaMin && aod::emcalcluster::eta < clusterEtaMax && aod::emcalcluster::phi >= clusterPhiMin && aod::emcalcluster::phi <= clusterPhiMax && aod::emcalcluster::energy >= clusterEnergyMin && aod::emcalcluster::time > clusterTimeMin && aod::emcalcluster::time < clusterTimeMax && (clusterRejectExotics && aod::emcalcluster::isExotic != true));


  // Process functions
  // Data tracks
  void processDataTracks(
    soa::Filtered<aod::JCollisions>::iterator const& collision,
    Tracks const& tracks)
  {
    if (!JetDerivedDataUtilities::selectCollision(collision, eventSelection)) {
      return;
    }

    // Event level

    //
    for (auto& track : tracks) {
      // NOTE: Charge is encoded into pt here!
      trackSkim(collision.globalIndex(), track.pt() * track.sign(), track.eta(), track.phi());
    }
  }
  PROCESS_SWITCH(TrackSkimExploratory, processDataTracks, "Data Tracks", true);

  // MC det and part level tracks
  void processMCTracks(
    soa::Join<aod::JCollisions, aod::JMcCollisionLbs>::iterator const& collision, soa::Join<aod::JTracks, aod::JMcTrackLabels> const& tracks,
    aod::JMcParticles const& mcParticles,
    aod::JMcCollisions const&)
  {
    // TODO: RJE: MC event selection
    // if (!JetDerivedDataUtilities::selectCollision(collision, eventSelection)) {
    //  return;
    //}

    // Event level
    if (!collision.has_mcCollision()) {
      return;
    }
    auto mcCollision = collision.mcCollision();
    eventSkimMC(collision.globalIndex(), mcCollision.weight());

    // Det level
    for (auto& track : tracks) {
      // FIXME: RJE: We will lose MC particles if unmatched. However, we don't have the right MC information with the match.
      //             This should be revisited later.
      if (!track.has_mcParticle()) {
        LOGF(warning, "No MC particle for track, skip...");
        continue;
      }
      auto particle = track.mcParticle();

      // NOTE: Charge is encoded into pt here!
      // TODO: RJE: Add encoded information when we decide what we want.
      trackSkimMCD(
        collision.globalIndex(),
        track.pt() * track.sign(), track.eta(), track.phi(),
        particle.pdgCode(), particle.globalIndex(), 0);
    }

    // Particle level
    for (auto& particle : mcParticles) {
      // NOTE: In principle, we don't need to encode the charge here since we have the PDG code.
      // TODO: RJE: Add encoded information when we decide what we want.
      trackSkimMCP(
        collision.globalIndex(),
        particle.pt(), particle.eta(), particle.phi(),
        particle.pdgCode(), particle.globalIndex(), 0);
    }
  }
  PROCESS_SWITCH(TrackSkimExploratory, processMCTracks, "Part and det level tracks", true);

  void processDataClusters(
    soa::Filtered<aod::JCollisions>::iterator const& collision,
    soa::Filtered<o2::aod::JCluster> const& clusters)
  {
    if (!JetDerivedDataUtilities::selectCollision(collision, eventSelection)) {
      return;
    }
    // Require EMCal
    if (!JetDerivedDataUtilities::eventEMCAL(collision)) {
      return;
    }

    for (auto& cluster : clusters) {
      clusterSkim(collision.globalIndex(), cluster.eta(), cluster.phi(), cluster.energy());
    }
  }
  PROCESS_SWITCH(TrackSkimExploratory, processDataClusters, "Data clusters", true);
  // TODO: Clusters at MC det level
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<TrackSkimExploratory>(cfgc, TaskName{"track-skim-exploratory"})};
}
