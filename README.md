# EleBremStudies

1) Install:

    * scram project CMSSW_13_X_Y
    * cd CMSSW_13_X_Y/src/
    * cmsenv
    * git cms-init
    * git clone git@github.com:bmarzocc/EleBremStudies.git
    * cd EleBremStudies
    * git checkout 13_X_Y
    * cd -
    * scram b -j 5

2) Run from miniAOD:

    * cd EleBremStudies/Dumpers/
    * cmsRun python/BremDumper_cfg.py

3) Run from RAW:

    * cd EleBremStudies/Dumpers/
    * cmsRun test/BremDumper_Data_MiniAODv2_GT36_fromRaw.py
