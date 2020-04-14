import os
import fnmatch
for i in os.listdir('/Volumes/MacOS_encrypted/Patient_data/In_process'):
        if os.path.isdir(os.path.join('/Volumes/MacOS_encrypted/Patient_data/In_process',i)):
            new_path = os.path.join('/Volumes/MacOS_encrypted/Patient_data/In_process',i)
            for j in os.listdir(new_path):
                if fnmatch.fnmatch(j,'SAGE.nii'):
                   file_name = os.path.join(new_path,j)
                  os.system('gzip %s' % file_name)
