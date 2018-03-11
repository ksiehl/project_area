# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: Hadronizer_TuneCUETP8M1_13TeV_generic_LHE_pythia8_cff.py --fileout file:output-aodsim.root --mc --eventcontent AODSIM --conditions auto:mc --step GEN,SIM,DIGI,L1,DIGI2RAW,HLT:@fake1,RAW2DIGI,RECO --python_filename config_pythia-aodsim.py --no_exec -n 1 --filein file:13TeVunweighted_events.lhe
import FWCore.ParameterSet.Config as cms

process = cms.Process('HLT')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.Geometry.GeometrySimDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedRealistic50ns13TeVCollision_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.Digi_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.DigiToRaw_cff')
process.load('HLTrigger.Configuration.HLT_Fake1_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(5)
)

# Input source
process.source = cms.Source("LHESource",
    dropDescendantsOfDroppedBranches = cms.untracked.bool(False),
    #fileNames = cms.untracked.vstring('file:13TeVunweighted_events.lhe'),
    fileNames = cms.untracked.vstring('file:dim4_0.lhe'),
    #fileNames = cms.untracked.vstring('file:sgnl-bckgrnd-stdy.lhe'),
    inputCommands = cms.untracked.vstring('keep *')#, 
        #'drop LHEXMLStringProduct_*_*_*', 
        #'drop *_genParticles_*_*', 
        #'drop *_genParticlesForJets_*_*', 
        #'drop *_kt4GenJets_*_*', 
        #'drop *_kt6GenJets_*_*', 
        #'drop *_iterativeCone5GenJets_*_*', 
        #'drop *_ak4GenJets_*_*', 
        #'drop *_ak7GenJets_*_*', 
        #'drop *_ak8GenJets_*_*', 
        #'drop *_ak4GenJetsNoNu_*_*', 
        #'drop *_ak8GenJetsNoNu_*_*', 
        #'drop *_genCandidatesForMET_*_*', 
        #'drop *_genParticlesForMETAllVisible_*_*', 
        #'drop *_genMetCalo_*_*', 
        #'drop *_genMetCaloAndNonPrompt_*_*', 
        #'drop *_genMetTrue_*_*', 
        #'drop *_genMetIC5GenJs_*_*')
)

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('Hadronizer_TuneCUETP8M1_13TeV_generic_LHE_pythia8_cff.py nevts:1'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.AODSIMoutput = cms.OutputModule("PoolOutputModule",
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('generation_step')
    ),
    compressionAlgorithm = cms.untracked.string('LZMA'),
    compressionLevel = cms.untracked.int32(4),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string(''),
        filterName = cms.untracked.string('')
    ),
    eventAutoFlushCompressedSize = cms.untracked.int32(15728640),
    fileName = cms.untracked.string('file:output/output-aodsim.root'),
    outputCommands = process.AODSIMEventContent.outputCommands
)

#####process.AODSIMoutput.outputCommands.append('keep *_allWeights_*_*') ##TEMPORARILY COMMENTING OUT WEIGHTING
#####process.AODSIMoutput.outputCommands.append('keep *_weightsMap_*_*') ##TEMPORARILY COMMENTING OUT WEIGHTING

# Additional output definition

# Other statements
process.genstepfilter.triggerConditions=cms.vstring("generation_step")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:mc', '')

process.generator = cms.EDFilter("Pythia8HadronizerFilter",
    PythiaParameters = cms.PSet(
        parameterSets = cms.vstring('pythia8CommonSettings', 
            'pythia8CUEP8M1Settings'),
        pythia8CUEP8M1Settings = cms.vstring('Tune:pp 14', 
            'Tune:ee 7', 
            'MultipartonInteractions:pT0Ref=2.4024', 
            'MultipartonInteractions:ecmPow=0.25208', 
            'MultipartonInteractions:expPow=1.6'),
        pythia8CommonSettings = cms.vstring('Tune:preferLHAPDF = 2', 
            'Main:timesAllowErrors = 10000', 
            'Check:epTolErr = 0.01', 
            'Beams:setProductionScalesFromLHEF = off', 
            'SLHA:keepSM = on', 
            'SLHA:minMassSM = 1000.', 
            'ParticleDecays:limitTau0 = on', 
            'ParticleDecays:tau0Max = 10', 
            'ParticleDecays:allowPhotonRadiation = on')
    ),
    comEnergy = cms.double(13000.0),
    filterEfficiency = cms.untracked.double(1.0),
    maxEventsToPrint = cms.untracked.int32(1),
    pythiaHepMCVerbosity = cms.untracked.bool(False),
    pythiaPylistVerbosity = cms.untracked.int32(1)
)

######adding multiple weights stuff
####from WSUAnalysis.LHEWeightProducer.lheWeightProducer_cfi import MGWeightsFromLHE##TEMPORARILY REMOVING WEIGHTS

####process.weightsMap = MGWeightsFromLHE.clone()##

####process.allWeights = MGWeightsFromLHE.clone(####
####    makeWeightsMap = cms.untracked.bool(False),
####    produceAllWeights = cms.untracked.bool(True),
####    )
####process.MGWeights = cms.Sequence(process.allWeights + process.weightsMap)####

######back to other code

process.ProductionFilterSequence = cms.Sequence(process.generator)

# Path and EndPath definitions
process.generation_step = cms.Path(process.pgen)
process.simulation_step = cms.Path(process.psim)
process.digitisation_step = cms.Path(process.pdigi)
process.L1simulation_step = cms.Path(process.SimL1Emulator)
process.digi2raw_step = cms.Path(process.DigiToRaw)# + process.MGWeights) ##TEMPORARILY REMOVING WEIGHTS
process.raw2digi_step = cms.Path(process.RawToDigi)
process.reconstruction_step = cms.Path(process.reconstruction)
process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.AODSIMoutput_step = cms.EndPath(process.AODSIMoutput)

# Schedule definition
process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step,process.simulation_step,process.digitisation_step,process.L1simulation_step,process.digi2raw_step)
process.schedule.extend(process.HLTSchedule)
process.schedule.extend([process.raw2digi_step,process.reconstruction_step,process.endjob_step,process.AODSIMoutput_step])
# filter all path with the production filter sequence
for path in process.paths:
	getattr(process,path)._seq = process.ProductionFilterSequence * getattr(process,path)._seq 

# customisation of the process.

# Automatic addition of the customisation function from HLTrigger.Configuration.customizeHLTforMC
from HLTrigger.Configuration.customizeHLTforMC import customizeHLTforFullSim 

#call to customisation function customizeHLTforFullSim imported from HLTrigger.Configuration.customizeHLTforMC
process = customizeHLTforFullSim(process)

# End of customisation functions

