# Medical_Image_Process_Fudan_University_School_of_Data_Science
复旦大学大数据学院2025春季医学图像处理
# description
project of 2025 Medical_Image_Process_Fudan_University_School_of_Data_Science
instructor is ZHOU Yuan
复旦大学大数据学院2025春季医学图像处理

# requirement
register spect to mri in mni space, using tpm.
using spm25

# procedure
## register each individual's spect to their own  mri
see coregister_job.m
## register each individual's mri to tpm
see normalize.m
## implement these two 
see deform job
## picture the output of the registration
included in reg.m
