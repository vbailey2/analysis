This is a package for looking at the fluctuations of the underlying event in Au+Au collisions using HIJING and eventually real data.

The code runs in two steps:
1. Tree making: this will run over dsts and produce trees containing kinematics for individual calorimeter towers in each event as well as information about the underlying event determined for the UE subtraction applied to jets. The code is based largely on the [CaloAna tutorial](https://github.com/sPHENIX-Collaboration/tutorials/tree/master/CaloAna). Instructions for building and running the code are below. 
  1. Build the package in the [usual way](https://wiki.bnl.gov/sPHENIX/index.php/Example_of_using_DST_nodes#Building%20a%20package):
   *Make a build directory inside the src directory: 
    cd src
    mkdir build
    cd build
   *Setup the sPHENIX environment and install paths:
    source /opt/sphenix/core/bin/sphenix_setup.csh -n
    setenv MYINSTALL ~/install
    source /opt/sphenix/core/bin/setup_local.csh $MYINSTALL
   *Compile your code and rerun the startup steps:
    make -j 4
    make install
    source /opt/sphenix/core/bin/sphenix_setup.csh -n
    setenv MYINSTALL ~/install
    source /opt/sphenix/core/bin/setup_local.csh $MYINSTALL
  2. Run the code using the Fun4All macro:
   *Go to the macro directory:
    cd ../../macro
   *Get some files to run on using CreateFileList.pl:
    CreateFileList.pl -type 4 -run 4 DST_CALO_CLUSTER
   *Test run using Fun4All. Note that the code is by default setup to read in files from dst_calo_cluster.list (the output from CreateFileList.pl and produce output file outputest.root, running over 1000 events. These file names can be specified as inputs to Fun4All_CaloAna() in case you want to specify them at run time, i.e. for running multiple jobs on condor:
    root -b -q Fun4All_CaloAna.C
   *This will create an output file containing all the necessary information for the histogram making in Step 2.
2. Histogram making: this will run over the output of Step 1 and produce histograms of the mean and standard deviation of windows of towers in each event vs. the event impact parameter and the total calorimeter energy. The method is described in slides [here](https://indico.bnl.gov/event/15444/contributions/62384/attachments/40602/67851/UEFluctuations_sPHENIXJetMeeting_4_20_22.pdf). To run:
  1. Go to the offline code directory- from the main directory (/analysis/UEFluctuations):
    cd offline
  2. Test run by doing:
    root -b -q getFluctuations.C
  3. Mess with the configuration. The default configuration will look at towers after the UE has been subtracted, however you can also study the unsubtracted towers, or towers that have the average UE subtracted without any corrections for flow. The inputs to getFluctuations() are as follows (string in = "../macro/outputest.root", string out ="testhist", int subtraction = 2, int nEvents = -1):
   *Input file (in): by default set to the file we produced in the test of ntuple making above
   *Output file (out): by default set to "testhist". Based on what you choose for the subtraction type, an additional string will be added to the end of the file name specifying subtracted/unsubtracted towers.
   *Subtraction type (subtraction): 0: no subtraction, 1: subtraction without flow, 2: subtraction with flow. By default set to 2.
   *Number of events (nEvents): for testing I'd suggest setting this to a small number, by default it is set to run over the full input file.
  4. Look at the output file: running the default configuration will produce the file testhist_subtracted_withflow.root. The titles of each histogram explain if they are the average or standard deviation, correspond to towers in a specific part of the detector (i.e. EMCal, inner or outer HCal) or the full calorimeter, if they are plotted vs. total calorimeter energy or impact parameter, and what window size they correspond to (1x1, 3x4, or 7x7).
  5. Draw nice plots. A script for doing this will be added soon.