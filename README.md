# fat-segmentation
In house MATLAB tool to semi-automatically segment fat/water-only DIXON scan images (some SPM8 dependency).

NOTE: to generate the rectified volumes from scratch requires some functions from SPM8.
This can be obtained from http://www.fil.ion.ucl.ac.uk/spm/software/download/.

Otherwise, simply use the subjectXX_Y_registered.mat provided as the files to load with VATsegmentation2.m.

The slice-by-slice segmentation results in subjectXX_Y_results.mat can be viewed with ShowResults.m. 
Comparison of pre- and post-intervention results can be viewed with CompareResutls.m.

-JEH
