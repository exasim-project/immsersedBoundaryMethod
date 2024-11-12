import sys, os 
import subprocess 

path=__file__.split('/compile')[0]

if not os.path.isdir(f'{path}/code/multiMultiFoam'):
    subprocess.run(['git', 'submodule', 'update', '--init'])
    subprocess.run(['git', 'apply', '../patch.multiMultiFoam'], cwd=f'{path}/code/multiMultiFoam')

print(f'{path}/code/multiMultiFoam')

subprocess.run(['./Allwmake'], cwd=f'{path}/code/multiMultiFoam')
subprocess.run(['wmake', '-j8', 'libso'], cwd=f'{path}/code/src')
