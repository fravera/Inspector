import FWCore.ParameterSet.Config as cms

process = cms.Process("DemoAnalysis")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
    noEventSort = cms.untracked.bool(True),
    fileNames = cms.untracked.vstring(
        #'file:/afs/cern.ch/work/f/fravera/Resonance/CMSSW_8_0_2/src/step_FastSimPPS_SIM_RECOBEFMIX_DIGI_L1_L1Reco_RECO_HLT_10kEvents.root'
        #'file:/afs/cern.ch/work/f/fravera/Resonance/NewFastSim/CMSSW_8_0_2/src/step_FastSimPPS_SIM_RECOBEFMIX_DIGI_L1_L1Reco_RECO_HLT_trackingOnly.root'
        # 'file:/afs/cern.ch/work/f/fravera/Resonance/NewFastSim/CMSSW_8_0_2/src/step_FastSimPPS_SIM_RECOBEFMIX_DIGI_L1_L1Reco_RECO_HLT_trackingWithTOF_10kEvents.root'
        'file:step_FastSimPPS_trackingWithToF_10kEvents.root'
        # 'file:/afs/cern.ch/work/f/fravera/Resonance/NewFastSim/CMSSW_8_0_2/src/step_FastSimPPS_SIM_RECOBEFMIX_DIGI_L1_L1Reco_RECO_HLT_prova_10Events.root'
        # 'root://xrootd.unl.edu//store/user/fravera/FPMC_aaaa_AAexoticSpin0EvenResonances_13TeV/CT_PPS_FastSim_Spin0Resonance/160511_164803/0000/step_FastSimPPS_SIM_RECOBEFMIX_DIGI_L1_L1Reco_RECO_HLT_1.root',
        # 'root://xrootd.unl.edu//store/user/fravera/FPMC_aaaa_AAexoticSpin0EvenResonances_13TeV/CT_PPS_FastSim_Spin0Resonance/160511_164803/0000/step_FastSimPPS_SIM_RECOBEFMIX_DIGI_L1_L1Reco_RECO_HLT_2.root',
        # 'root://xrootd.unl.edu//store/user/fravera/FPMC_aaaa_AAexoticSpin0EvenResonances_13TeV/CT_PPS_FastSim_Spin0Resonance/160511_164803/0000/step_FastSimPPS_SIM_RECOBEFMIX_DIGI_L1_L1Reco_RECO_HLT_3.root',
        # 'root://xrootd.unl.edu//store/user/fravera/FPMC_aaaa_AAexoticSpin0EvenResonances_13TeV/CT_PPS_FastSim_Spin0Resonance/160511_164803/0000/step_FastSimPPS_SIM_RECOBEFMIX_DIGI_L1_L1Reco_RECO_HLT_4.root',
        # 'root://xrootd.unl.edu//store/user/fravera/FPMC_aaaa_AAexoticSpin0EvenResonances_13TeV/CT_PPS_FastSim_Spin0Resonance/160511_164803/0000/step_FastSimPPS_SIM_RECOBEFMIX_DIGI_L1_L1Reco_RECO_HLT_5.root',
        # 'root://xrootd.unl.edu//store/user/fravera/FPMC_aaaa_AAexoticSpin0EvenResonances_13TeV/CT_PPS_FastSim_Spin0Resonance/160511_164803/0000/step_FastSimPPS_SIM_RECOBEFMIX_DIGI_L1_L1Reco_RECO_HLT_6.root',
        # 'root://xrootd.unl.edu//store/user/fravera/FPMC_aaaa_AAexoticSpin0EvenResonances_13TeV/CT_PPS_FastSim_Spin0Resonance/160511_164803/0000/step_FastSimPPS_SIM_RECOBEFMIX_DIGI_L1_L1Reco_RECO_HLT_7.root',
        # 'root://xrootd.unl.edu//store/user/fravera/FPMC_aaaa_AAexoticSpin0EvenResonances_13TeV/CT_PPS_FastSim_Spin0Resonance/160511_164803/0000/step_FastSimPPS_SIM_RECOBEFMIX_DIGI_L1_L1Reco_RECO_HLT_8.root',
        # 'root://xrootd.unl.edu//store/user/fravera/FPMC_aaaa_AAexoticSpin0EvenResonances_13TeV/CT_PPS_FastSim_Spin0Resonance/160511_164803/0000/step_FastSimPPS_SIM_RECOBEFMIX_DIGI_L1_L1Reco_RECO_HLT_9.root',
        # 'root://xrootd.unl.edu//store/user/fravera/FPMC_aaaa_AAexoticSpin0EvenResonances_13TeV/CT_PPS_FastSim_Spin0Resonance/160511_164803/0000/step_FastSimPPS_SIM_RECOBEFMIX_DIGI_L1_L1Reco_RECO_HLT_10.root'
    )
)

process.demo = cms.EDAnalyzer('InspectReco'
)

process.TFileService = cms.Service("TFileService",
                                       #fileName = cms.string('/afs/cern.ch/work/f/fravera/Resonance/NewFastSim/CMSSW_8_0_2/src/histoInspectReco_trackingOnly.root')
                                       fileName = cms.string('histoInspectReco_trackingWithToF_10kEvents.root')
                                       # fileName = cms.string('/afs/cern.ch/work/f/fravera/Resonance/NewFastSim/CMSSW_8_0_2/src/histoInspectReco_prova_10Events.root')
                                       #fileName = cms.string('/afs/cern.ch/work/f/fravera/Resonance/CMSSW_8_0_2/src/histoInspectReco_trackingOnly.root')
                                   )

process.p = cms.Path(process.demo)
