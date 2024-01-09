

# CRAB3 config template for flashgg
# More options available on the twiki :
# https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookCRAB3Tutorial

from WMCore.Configuration import Configuration

config = Configuration()

config.section_('General')
config.General.requestName       = 'EGamma_Run2018A-UL2018_MiniAODv2_GT36-v1'
config.General.transferLogs      = True
config.General.transferOutputs   = True

config.section_('JobType')
config.JobType.pluginName        = 'Analysis'

# Name of the CMSSW configuration file
config.JobType.psetName          = 'BremDumper_Data_MiniAODv2_GT36_fromRaw.py'
#config.JobType.psetName          = 'BremDumper_cfg.py'
config.JobType.priority          = 30
config.JobType.maxMemoryMB       = 5500
config.JobType.numCores          = 4
config.JobType.inputFiles        = ['Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt'] 
config.JobType.outputFiles       = ['output.root']
config.JobType.allowUndistributedCMSSW = True
config.JobType.disableAutomaticOutputCollection = True

config.section_('Data')
# This string determines the primary dataset of the newly-produced outputs.

config.Data.inputDataset          = '/EGamma/Run2018A-UL2018_MiniAODv2_GT36-v1/MINIAOD'
config.Data.secondaryInputDataset = '/EGamma/Run2018A-v1/RAW'

#config.Data.inputDataset          = '/EGamma/Run2018B-UL2018_MiniAODv2_GT36-v1/MINIAOD'
#config.Data.secondaryInputDataset = '/EGamma/Run2018B-v1/RAW'

#config.Data.inputDataset          = '/EGamma/Run2018C-UL2018_MiniAODv2_GT36-v1/MINIAOD'
#config.Data.secondaryInputDataset = '/EGamma/Run2018C-v1/RAW'

#config.Data.inputDataset          = '/EGamma/Run2018D-UL2018_MiniAODv2_GT36-v3/MINIAOD'
#config.Data.secondaryInputDataset = '/EGamma/Run2018D-v1/RAW'

config.Data.inputDBS             = 'global'   
#config.Data.inputDBS             = 'phys03'   
config.Data.splitting            = 'FileBased'
config.Data.unitsPerJob          = 1
config.Data.lumiMask             = 'Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt'
config.Data.publication          = False
#config.Data.ignoreLocality       = True
config.Data.outputDatasetTag     = 'Run2018A-UL2018'

# This string is used to construct the output dataset name
#config.Data.publishDataName = 'CRAB3-tutorial'
config.Data.outLFNDirBase        =  '/store/group/phys_higgs/cmshgg/bmarzocc/FNUF/PhotonEoPDumper'

config.section_('Site')
# Where the output files will be transmitted to
config.Site.storageSite         = 'T2_CH_CERN'
#config.Site.whitelist           = ['T2_CH_CERN']


## config.Data.allowNonValidInputDataset=True
