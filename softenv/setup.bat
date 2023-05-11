REM set up conda environments needed for processing barseq

mkdir c:\barseq_envs\
echo Extracting environments from archive...
tar -xf bardensr.tar -C C:\barseq_envs\
tar -xf cellpose.tar -C C:\barseq_envs\
tar -xf n2v.tar -C C:\barseq_envs\
rem append additional tars here if needed



rem add the environment folder to envs_dirs
echo Adding environment path to conda envs_dirs
conda config --append envs_dirs c:\barseq_envs
echo Done!
