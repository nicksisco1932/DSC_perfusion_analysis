import fnmatch
import multiprocessing as mp
import os
import shutil
import dicom2nifti
from tkinter import filedialog
from tkinter import *

T1w_pre_dir_names = ['WIP AX 3D T1', 'WIP AX T1 3D WAND', '*AX T1 3D_TFE (WAND)*']
T1w_post_dir_names = ['WIP +CsT1W_3D_TFE*', 'WIP AX T1 3D TFE+C', '*+CsT1W_3D_TFE*', '*AX T1 3D_TFE (WAND)+C*',
                      'WIP AX T1 3D_TFE (WAND)+C 1_2']
# T2w_FLAIR_dir_names = ['*AX*FLAIR*', 'WIP AX T2 FLAIR MVXD', 'WIP AX FLAIR MVXD']
T2w_FLAIR_dir_names = ['*AX*FLAIR*']
SAGE_dir_names = ['WIP INJECT_SAGEperfusion*', '*INJECT_SAGEperfusion*', 'SAGE_perfusion']

in_names = ["T1w_pre.nii", "T1w_post.nii", "T2w_FLAIR.nii", "SAGE.nii"]
out_names = ["brainT1w_pre.nii", "brainT1w_post.nii", "brainT2w_FLAIR.nii", "brainSAGE.nii"]
zip_files = ['brainT1w_pre.nii.gz', 'brainT1w_post.nii.gz', 'brainT2w_FLAIR.nii.gz', 'brainSAGE.nii.gz']


def main():
    root = Tk()
    root.withdraw()
    path = filedialog.askdirectory()
    print(path)
    work(path)


def bet_proc(path, filename):
    # in_names = ["T1w_pre.nii", "T1w_post.nii", "T2w_FLAIR.nii", "SAGE.nii"]
    # out_names = ["brainT1w_pre.nii", "brainT1w_post.nii", "brainT2w_FLAIR.nii", "brainSAGE.nii"]
    # zip_files = ['brainT1w_pre.nii.gz', 'brainT1w_post.nii.gz', 'brainT2w_FLAIR.nii.gz', 'brainSAGE.nii.gz']
    temp = path
    run_once = 0
    if os.path.isdir(path):
        for j in [0, 1]:
            if not os.path.isfile(os.path.join(path, in_names[j])):
                print(in_names[j] + ' does not exist')
                pass
            else:
                if not os.path.isfile(os.path.join(path, out_names[j])) and not os.path.isfile(
                        os.path.join(path, zip_files[j])):
                    # os.system("bet %s %s -f 0.5 -g -0.4 -m > /dev/null" % (in_names[j], out_names[j]))
                    os.system('bet %s %s -f 0.5 -g -0.4 -m ' % (
                        os.path.join(path, in_names[j]), os.path.join(path, out_names[j])))
                    # os.system("gunzip -f %s" % zip_files[j])
                else:
                    pass
        for j in [2, 3]:
            if not os.path.isfile(os.path.join(path, in_names[j])):
                print(in_names[j] + ' does not exist')
                pass
            else:
                if not os.path.isfile(os.path.join(path, out_names[j])) and not os.path.isfile(
                        os.path.join(path, zip_files[j])):
                    # os.system("bet %s %s -f 0.4 -g 0 -m > /dev/null" % (in_names[j], out_names[j]))
                    os.system('bet %s %s -f 0.4 -g -0 -m ' % (
                        os.path.join(path, in_names[j]), os.path.join(path, out_names[j])))
                    # os.system("gunzip -f %s" % zip_files[j])
                else:
                    pass
    else:
        pass
    print('Brain Extraction and Brain Mask Creation Finished for ' + filename)


def work(path):
    pool1 = mp.Pool(processes=os.cpu_count())
    for i in os.listdir(path):
        if os.path.isdir(os.path.join(path, i)):
            path_1 = os.path.join(path, i)
            # print(path_1)
            print(i)
            p1 = pool1.apply_async(bet_proc, args=(path_1, i,))
            p2 = pool1.apply_async(register_everything, args=(path_1, i,))
            # bet_proc(path_1, i)
            # register_everything(path_1, i)
        else:
            pass
    pool1.close()  # do not accept any more tasks
    pool1.join()


def register_everything(path, filename):
    # This does all of the registration. The first loop is registering T1wPost to T1wPre, the second loop registers
    # T1wPost and T1wPre to the registration of the SAGE experiments, and the last loop registers T2wFLAIR to SAGE
    # # resolution.
    # in_names = [os.path.join(path, "T1w_pre.nii"), os.path.join(path, "T1w_post.nii"),
    #             os.path.join(path, "T2w_FLAIR.nii"),
    #             os.path.join(path, "SAGE.nii")]
    # out_names = [os.path.join(path, "brainT1w_pre.nii"), os.path.join(path, "brainT1w_post.nii"),
    #              os.path.join(path, "brainT2w_FLAIR.nii"), os.path.join(path, "brainSAGE.nii")]
    # zip_files = [os.path.join(path, 'brainT1w_pre.nii.gz'), os.path.join(path, 'brainT1w_post.nii.gz'),
    #              os.path.join(path, 'brainT2w_FLAIR.nii.gz'), os.path.join(path, 'brainSAGE.nii.gz')]

    if not os.path.isfile(os.path.join(path, 'AffineReg12p_brainT1wpost2pre.nii.gz')):
        if not os.path.isfile(os.path.join(path, out_names[1] + '.gz')):
            pass
        else:
            os.system('flirt -bins 256 -cost normmi -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -dof 12 '
                      '-interp trilinear -in %s -out %s'
                      ' -ref %s -omat %s' % (os.path.join(path, out_names[1] + '.gz'),
                                             os.path.join(path, 'AffineReg12p_brainT1wpost2pre.nii'),
                                             os.path.join(path, out_names[0] + '.gz'),
                                             os.path.join(path, 'Reg_refT1wpre_inputT1wpost.txt')))
    elif not os.path.isfile(os.path.join(path, 'AffineReg12p_brainT1wpre2pre_SAGEres.nii.gz')):
        if not os.path.isfile(os.path.join(path, out_names[1] + '.gz')):
            pass
        else:
            os.system('flirt -bins 256 -cost normmi -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -dof 12 '
                      '-interp trilinear -in %s -out %s'
                      ' -ref %s -omat %s' % (os.path.join(path, out_names[0] + '.gz'),
                                             os.path.join(path, 'AffineReg12p_brainT1wpre2pre_SAGEres.nii'),
                                             os.path.join(path, out_names[3] + '.gz'),
                                             os.path.join(path, 'Reg_refT1wpreSAGE_inputT1wpre.txt')))
    elif not os.path.isfile(os.path.join(path, 'T1w_pre_SAGEres.nii.gz')):
        if not os.path.isfile(in_names[0]):
            pass
        else:
            os.system('flirt -paddingsize 0.0 -interp trilinear -in %s -out %s '
                      '-ref %s -applyxfm -init %s' % (os.path.join(path, in_names[0]),
                                                      os.path.join(path, 'T1w_pre_SAGEres.nii'),
                                                      os.path.join(path, out_names[3] + '.gz'),
                                                      os.path.join(path, 'Reg_refT1wpreSAGE_inputT1wpre.txt')))
    elif not os.path.isfile(os.path.join(path, 'brainT1w_pre_SAGEres.nii.gz')):
        if not os.path.isfile(os.path.join(path, out_names[0] + '.gz')):
            pass
        else:
            os.system('flirt -paddingsize 0.0 -interp trilinear -in %s -out %s '
                      '-ref %s' % (os.path.join(path, out_names[0] + '.gz'),
                                   os.path.join(path, 'brainT1w_pre_SAGEres'),
                                   os.path.join(path, out_names[3] + '.gz'))),

    elif not os.path.isfile(os.path.join(path, 'AffineReg12p_brainT1wpost2pre_SAGEres.nii.gz')):
        if not os.path.isfile(os.path.join(path, out_names[1] + '.gz')):
            pass
        else:
            os.system('flirt -bins 256 -cost normmi -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -dof 12 '
                      '-interp trilinear -in %s -out %s'
                      ' -ref %s -omat %s'
                      % (os.path.join(path, out_names[1] + '.gz'),
                         os.path.join(path, 'AffineReg12p_brainT1wpost2pre_SAGEres'),
                         os.path.join(path, 'brainT1w_pre_SAGEres.nii.gz'),
                         os.path.join(path, 'Reg_refT1wpreSAGE_inputT1wpost.txt')))

    elif not os.path.isfile(os.path.join(path, 'T1w_post_SAGEres.nii.gz')):
        if not os.path.isfile(os.path.join(path, in_names[1])):
            pass
        else:
            os.system('flirt -paddingsize 0.0 -interp trilinear -in %s -out %s '
                      '-ref %s -applyxfm -init %s' % (os.path.join(path, in_names[1]),
                                                      os.path.join(path, 'T1w_post_SAGEres'),
                                                      os.path.join(path, out_names[3] + '.gz'),
                                                      os.path.join(path, 'Reg_refT1wpreSAGE_inputT1wpost.txt')))
    elif not os.path.isfile(os.path.join(path, 'AffineReg12p_brainSAGE.nii.gz')):
        if not os.path.isfile(os.path.join(path, out_names[3] + '.gz')):
            pass
        else:
            os.system('flirt -bins 256 -cost normmi -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -dof 12 '
                      '-interp trilinear -in %s -out %s'
                      ' -ref %s -omat %s' % (os.path.join(path, out_names[3] + '.gz'),
                                             os.path.join(path, 'AffineReg12p_brainSAGE'),
                                             os.path.join(path, out_names[2] + '.gz'),
                                             os.path.join(path, 'Reg_refT2w_inputSAGE.txt')))
            os.system('convert_xfm -omat %s -inverse %s'
                      % (os.path.join(path, 'Reg_refSAGE_inputT2w_FLAIR.txt'),
                         os.path.join(path, 'Reg_refT2w_inputSAGE.txt')))

    elif not os.path.isfile(os.path.join(path, 'T2w_FLAIR_SAGEres.nii.gz')):
        if not os.path.isfile(os.path.join(path, in_names[2])):
            pass
        else:
            os.system('flirt -paddingsize 0.0 -interp trilinear -in %s -out %s '
                      '-ref %s -applyxfm -init %s' % (os.path.join(path, in_names[2]),
                                                      os.path.join(path, 'T2w_FLAIR_SAGEres'),
                                                      os.path.join(path, out_names[3] + '.gz'),
                                                      os.path.join(path, 'Reg_refSAGE_inputT2w_FLAIR.txt')))
    elif not os.path.isfile(os.path.join(path, 'brainT2w_FLAIR_SAGEres.nii.gz')):
        if not os.path.isfile(os.path.join(path, out_names[2] + '.gz')):
            pass
        else:
            os.system('flirt -paddingsize 0.0 -interp trilinear -in %s -out %s '
                      '-ref %s -applyxfm -init %s' % (os.path.join(path, out_names[2] + '.gz'),
                                                      os.path.join(path, 'brainT2w_FLAIR_SAGEres.nii.gz'),
                                                      os.path.join(path, out_names[3]) + '.gz',
                                                      os.path.join(path, 'Reg_refSAGE_inputT2w_FLAIR.txt')))
    else:
        pass
    print('Registration Finished for ' + filename)


if __name__ == '__main__':
    main()
