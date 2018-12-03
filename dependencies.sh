# Set fsl version to 5.0.11
set_fsl_version 5.0.11

# Set conda environment to the one that has neuroprocess functions
export PATH=/study/midus3/environments/conda3/bin:$PATH
export PYTHONHOME=/study/midus3/environments/conda3

# Set path to BSE
LD_LIBRARY_PATH=/apps/MCR-R2015b/x86_64_sci7/v90/runtime/glnxa64:/apps/MCR-R2015b/x86_64_sci7/v90/bin/glnxa64:/apps/MCR-R2015b/x86_64_sci7/v90/sys/os/glnxa64
PATH=/apps/BrainSuite-18a/x86_64_sci7/bin:$PATH

export PS1="(midus3-iamge-processing) $PS1"