from WMCore.Configuration import Configuration
config = Configuration()
config.section_('General')
config.General.transferOutputs = True
config.section_('JobType')
config.JobType.psetName = 'config_pythia-aodsim.py'
config.JobType.pluginName = 'privateMC'
#config.JobType.inputFiles = ['/uscms_data/d3/ksiehl/cp_violating_anom/LHE_root_conversion_source_files/sgnl-bckgrnd-stdy.lhe']
config.JobType.inputFiles = ['dim4_0.lhe']# ['13TeVunweighted_events.lhe']
config.JobType.generator = 'lhe'
config.section_('Data')
config.Data.outputPrimaryDataset = 'test_run'
#config.Data.inputDataset = 'None'
#config.Data.publishDBS = 'https://cmsdbsprod.cern.ch:8443/cms_dbs_ph_analysis_02_writer/servlet/DBSServlet'
config.Data.publication = True
config.Data.unitsPerJob = 1500
config.Data.splitting = 'EventBased'
config.Data.outputDatasetTag = 'test-run'
config.Data.totalUnits = 50000
config.section_('User')
config.section_('Site')
config.Site.storageSite = 'T3_US_FNALLPC'
