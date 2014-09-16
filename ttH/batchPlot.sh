qsub -N ttH  -o /afs/desy.de/user/e/eckardt/Plotting/ttH/logs/ -e /afs/desy.de/user/e/eckardt/Plotting/ttH/logs/ -m eas -l cvmfs -l h_rt=24:00:00 -l h_vmem=12G -l h_fsize=12G batchPlotScript.sh
