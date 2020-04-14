# DSC_perfusion_analysis

% The whole idea behind this project here is to make the data uniform across the patients so the analysis is easier later on. 

% Created by Nicholas J. Sisco, Ph.D. please send questions to
% nicholas dot sisco at barrowneuro dot org

One major caveat to this project is that the data has to be sorted using eith Horos software or https://dicomsort.com.

It's a work in progress but the general flow goes:
1) Perfusion_analysis.m
2) SAGE_proc_only.m
3) DSC_registration_NS.py

Later versions will have SAGE_proc_only.m nested in Perfusion_analysis.m

The  DSC_registration_NS.py is really useful as it is multiprocessor Python script that makes light duty of registration using FSL.

NOTE: You'll need to gzip all of the SAGE.nii files for the Python script to work. Use this. 
for i in os.listdir('/Volumes/MacOS_encrypted/Patient_data/In_process'):
        if os.path.isdir(os.path.join('/Volumes/MacOS_encrypted/Patient_data/In_process',i)):
            new_path = os.path.join('/Volumes/MacOS_encrypted/Patient_data/In_process',i)
            for j in os.listdir(new_path):
                if fnmatch.fnmatch(j,'SAGE.nii'):
                   file_name = os.path.join(new_path,j)
                  os.system('gzip %s' % file_name)

