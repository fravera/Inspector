import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        # 'root://xrootd.unl.edu//store/user/antoniov/FPMC_aaaa_AAexoticSpin0EvenResonances_13TeV/FPMC_aaaa_AAexoticSpin0EvenResonances_13TeV-v2/160504_150338/0000/FPMC_HepMC_GEN_1.root'
        'file:FPMC_aaaa_spin0_13TeV_HepMC_GEN.root'
        #'file:4875BCA4-E5CF-E311-B217-003048D3750A.root'
        #'root://eoscms.cern.ch//eos/cms/store/mc/RunIIFall15MiniAODv2/GluGluHToGG_M-125_13TeV_powheg_pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/50000/02ABA756-B9B8-E511-9578-002618943C3A.root'
        #'/store/generator/Summer13/HToGammaGamma_125_14TeV_powheg_pythia6/GEN/UpgrdPhase2BE_POSTLS261_V2-v1/00000/56677DFA-8CC7-E211-A7A6-02163E008C0F.root'
        #'root://xrootd.unl.edu//store/mc/RunIISummer15GS/THW_HToGG_13TeV-madgraph-pythia8_TuneCUETP8M1/GEN-SIM/MCRUN2_71_V1-v2/10000/00973E66-FAEB-E511-B9F4-02163E00F304.root'
        #'root://xrootd.unl.edu//store/mc/RunIISummer15GS/THW_HToGG_13TeV-madgraph-pythia8_TuneCUETP8M1/GEN-SIM/MCRUN2_71_V1-v2/10000/2E4B37A4-DCEB-E511-918B-02163E00E7FE.root'
    )
)

process.demo = cms.EDAnalyzer('InspectGenerator'
)

process.TFileService = cms.Service("TFileService",
                                       fileName = cms.string('histoInspectGenerator.root')
                                   )

process.p = cms.Path(process.demo)
