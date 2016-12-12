import FWCore.ParameterSet.Config as cms

process = cms.Process("Ntupler")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10

process.load("Geometry.CMSCommonData.cmsIdealGeometryXML_cfi");
process.load("Geometry.CaloEventSetup.CaloGeometry_cfi");
process.load("Geometry.CaloEventSetup.CaloTopology_cfi");

#
# Define input data to read
#
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100000) )

#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#process.GlobalTag.globaltag = 'MCRUN2_74_V9::All'                                             #MC
#process.GlobalTag.globaltag = 'GR_E_V49::All'                                                 #Data

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
#process.GlobalTag.globaltag = cms.string("76X_mcRun2_asymptotic_v12")
process.GlobalTag.globaltag = cms.string("76X_dataRun2_v15")
process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
        calibratedPatElectrons = cms.PSet(
        initialSeed = cms.untracked.uint32(1),
        engineName = cms.untracked.string('TRandom3')
        ),
        calibratedElectrons = cms.PSet(
        initialSeed = cms.untracked.uint32(1),
        engineName = cms.untracked.string('TRandom3')
        ),
                                                   )

inputFilesAOD = cms.untracked.vstring(
    # AOD test files from /DYJetsToLL_M-50_13TeV-madgraph-pythia8/Phys14DR-PU20bx25_PHYS14_25_V1-v1/AODSIM
    '/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/00CC714A-F86B-E411-B99A-0025904B5FB8.root',
    )    

inputFilesMiniAOD = cms.untracked.vstring(
    # MiniAOD test files from /DYJetsToLL_M-50_13TeV-madgraph-pythia8/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIM

#'/store/data/Run2015D/SingleElectron/MINIAOD/16Dec2015-v1/20000/00050EF1-F9A6-E511-86B2-0025905A48D0.root'
#'/store/mc/RunIIFall15MiniAODv2/ZToEE_NNPDF30_13TeV-powheg_M_2300_3500/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/60000/C84500E1-6BB8-E511-A9D6-002590E505FE.root'

       '/store/mc/RunIIFall15MiniAODv2/DYJetsToLL_M-200to400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext1-v1/10000/02DC197B-9AC9-E511-A223-0025907D244A.root',
       '/store/mc/RunIIFall15MiniAODv2/DYJetsToLL_M-200to400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext1-v1/10000/34D22204-9AC9-E511-908B-00259029ED16.root',
       '/store/mc/RunIIFall15MiniAODv2/DYJetsToLL_M-200to400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext1-v1/10000/4A556A46-9AC9-E511-B255-0CC47A0AD6EA.root',
       '/store/mc/RunIIFall15MiniAODv2/DYJetsToLL_M-200to400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext1-v1/10000/648B2234-9AC9-E511-9567-00259029E720.root',
       '/store/mc/RunIIFall15MiniAODv2/DYJetsToLL_M-200to400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext1-v1/10000/8A92AB05-9AC9-E511-A82D-00259048A862.root',
       '/store/mc/RunIIFall15MiniAODv2/DYJetsToLL_M-200to400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext1-v1/10000/9447F704-9AC9-E511-8BC1-0025904CBF10.root',
       '/store/mc/RunIIFall15MiniAODv2/DYJetsToLL_M-200to400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext1-v1/10000/B2D65535-9AC9-E511-AF1E-00259029E84C.root',
       '/store/mc/RunIIFall15MiniAODv2/DYJetsToLL_M-200to400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext1-v1/10000/D2CB6B46-9AC9-E511-9F7B-0CC47A0AD6EA.root',
       '/store/mc/RunIIFall15MiniAODv2/DYJetsToLL_M-200to400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext1-v1/10000/FA5C6629-9AC9-E511-B8FA-0CC47A0AD6EA.root',
       '/store/mc/RunIIFall15MiniAODv2/DYJetsToLL_M-200to400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext1-v1/60000/1653D32B-1FC8-E511-88F4-002590D9D9DA.root',
       '/store/mc/RunIIFall15MiniAODv2/DYJetsToLL_M-200to400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext1-v1/60000/28346BD3-ABC7-E511-8D27-001E67DFF6E0.root',
       '/store/mc/RunIIFall15MiniAODv2/DYJetsToLL_M-200to400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext1-v1/60000/3A1831CC-19C8-E511-B75F-00304867FE73.root',
       '/store/mc/RunIIFall15MiniAODv2/DYJetsToLL_M-200to400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext1-v1/60000/4204342F-1FC8-E511-B613-003048CBA444.root',
       '/store/mc/RunIIFall15MiniAODv2/DYJetsToLL_M-200to400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext1-v1/60000/76E9BC2D-24C8-E511-8DE4-0025901ABD1A.root',
       '/store/mc/RunIIFall15MiniAODv2/DYJetsToLL_M-200to400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext1-v1/60000/928D53F2-E5C7-E511-ABF9-002264031B48.root'

)
useAOD = False
if useAOD == True :
    inputFiles = inputFilesAOD
    print("AOD input files are used")
else :
    inputFiles = inputFilesMiniAOD
    print("MiniAOD input files are used")
process.source = cms.Source ("PoolSource", fileNames = inputFiles )                             



#process.load('EgammaAnalysis.ElectronTools.calibratedElectronsRun2_cfi')


# Set up electron ID (VID framework)

from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
# turn on VID producer, indicate data format  to be
# DataFormat.AOD or DataFormat.MiniAOD, as appropriate 
if useAOD == True :
    dataFormat = DataFormat.AOD
else :
    dataFormat = DataFormat.MiniAOD

switchOnVIDElectronIdProducer(process, dataFormat)

# define which IDs we want to produce
my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Spring15_25ns_V1_cff',
#                 'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring15_25ns_Trig_V1_cff',
#	         'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring15_25ns_nonTrig_V1_cff',
                 'RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV60_cff']


#add them to the VID producer
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)

#process.egmGsfElectronIDs.physicsObjectSrc = cms.InputTag('calibratedPatElectrons')
#process.egmGsfElectronIDs.physicsObjectSrc = cms.InputTag('selectedElectrons')
#
# Configure the ntupler module
#

process.selectedElectrons = cms.EDFilter("PATElectronSelector", 
src = cms.InputTag("slimmedElectrons"), 
cut = cms.string("pt > 10 && abs(eta)<2.5") 
)

process.ntupler = cms.EDAnalyzer('SimpleElectronNtupler',
                                 # The module automatically detects AOD vs miniAOD, so we configure both
                                 #
                                 # Common to all formats objects
                                 #
                                 isMC     = cms.untracked.bool(True),
				 isSIG    = cms.untracked.bool(True),
				 trigger  = cms.InputTag("TriggerResults::HLT"),
				 prescale = cms.InputTag("patTrigger"),
				 #pileup   = cms.InputTag("addPileupInfo"),
                                 pileup   = cms.InputTag("slimmedAddPileupInfo"),
				 rho      = cms.InputTag("fixedGridRhoFastjetAll"),
                                 beamSpot = cms.InputTag('offlineBeamSpot'),
                                 eventWeight   = cms.InputTag("generator"),
                                 #
                                 # Objects specific to AOD format
                                 #
 #                                electrons    = cms.InputTag("gedGsfElectrons"),
                                 electrons    = cms.InputTag("calibratedElectrons"),
                                 genParticles = cms.InputTag("genParticles"),
                                 vertices     = cms.InputTag("offlinePrimaryVertices"),
                                 conversions  = cms.InputTag('allConversions'),
                                 #
                                 # Objects specific to MiniAOD format
                                 #
                                 muonsMiniAOD    = cms.InputTag("slimmedMuons"),
				 electronsMiniAOD    = cms.InputTag("slimmedElectrons"),
    #                            electronsMiniAOD    = cms.InputTag("selectedElectrons"),
#                                 electronsMiniAOD    = cms.InputTag("calibratedPatElectrons"),
				 photonsMiniAOD      = cms.InputTag("slimmedPhotons"),
				 metsMiniAOD         = cms.InputTag("slimmedMETs"),
                                 genParticlesMiniAOD = cms.InputTag("prunedGenParticles"),
                                 verticesMiniAOD     = cms.InputTag("offlineSlimmedPrimaryVertices"),
                                 conversionsMiniAOD  = cms.InputTag('reducedEgamma:reducedConversions'),
				 #
                                 # ID decisions (common to all formats)
                                 #
				 eleVetoIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-veto"),
                                 eleLooseIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-loose"),
                                 eleMediumIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-medium"),
                                 eleTightIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-tight"),
				 eleHEEPIdMap = cms.InputTag("egmGsfElectronIDs:heepElectronID-HEEPV60"),
 #                                eletrigMVAlooseIdMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring15-25ns-Trig-V1-wp90"),
 ##                                eletrigMVAtightIdMap =  cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring15-25ns-Trig-V1-wp80"),
  #                               elenontrigMVAlooseIdMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring15-25ns-nonTrig-V1-wp90"),
  #                               elenontrigMVAtightIdMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring15-25ns-nonTrig-V1-wp80"),
				 phoLooseIdMap = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Spring15-25ns-V1-standalone-loose"),
                                 phoMediumIdMap = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Spring15-25ns-V1-standalone-medium"),
                                 phoTightIdMap = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Spring15-25ns-V1-standalone-tight"),
				 #
                                 # ValueMaps with MVA results
                                 #
                                 #mvaValuesMap     = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Phys14NonTrigValues"),
                                 #mvaCategoriesMap = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Phys14NonTrigCategories"),
				 # This is a fairly verbose mode if switched on, with full cut flow 
				 # diagnostics for each candidate. Use it in a low event count test job.
				 eleMediumIdFullInfoMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-medium",),
				 eleIdVerbose = cms.bool(False),
				 objects = cms.InputTag('selectedPatTrigger')
				 #tagFilterName   = cms.vstring("hltEle23Ele12CaloIdTrackIdIsoTrackIsoLeg1Filter"),
				 #probeFilterName = cms.vstring("hltEle23Ele12CaloIdTrackIdIsoTrackIsoLeg2Filter")
                                 )

process.primaryVertexFilter  = cms.EDFilter("VertexSelector",
      src = cms.InputTag('offlineSlimmedPrimaryVertices'),
      cut = cms.string('!isFake && ndof > 4.0 && position.Rho < 2.0 && abs(z) < 24'),
      filter = cms.bool(True)  ## otherwise it won't filter the events, just produce an empty vertex collection.
      )

process.hcalDDDRecConstants = cms.ESProducer( "HcalDDDRecConstantsESModule",
  appendToDataLabel = cms.string( "" )
)
process.hcalDDDSimConstants = cms.ESProducer( "HcalDDDSimConstantsESModule",
  appendToDataLabel = cms.string( "" )
)

#from PhysicsTools.PatAlgos.preselection_eltau_cff import *
#process.leadingElectron = process.slimmedElectrons.clone(
#    src = 'slimmedElectrons', cut = 'pt>20'
#    )
#process.load("PhysicsTools.PatAlgos.selectionLayer1.electronCountFilter_cfi")
#process.leadingElectronRequirement = process.countPatElectrons.clone(minNumber = 1, src = 'leadingElectron')
#process.leadingleptonsequence = cms.Sequence(process.leadingElectron+process.leadingElectronRequirement)

process.TFileService = cms.Service("TFileService",
 #                                  fileName = cms.string('DY_MCNLO_50toInf_76X.root')
				   fileName = cms.string('Data.root')
                                   )


#process.p = cms.Path(process.selectedElectrons + process.calibratedPatElectrons + process.egmGsfElectronIDSequence + process.primaryVertexFilter + process.ntupler) 
process.p = cms.Path(process.egmGsfElectronIDSequence * process.primaryVertexFilter * process.ntupler)
