import subprocess
import os, sys
import numpy as np

venvName='imageVenv'

if not os.path.isdir(venvName):
    print(f'Creating venv using {sys.executable}')    
    subprocess.run([sys.executable, '-m', 'venv', venvName], check=True)
    subprocess.run([f'{venvName}/bin/python', '-m', 'pip', 'install', '-r', 'requirements.txt'], check=True)

if not os.path.isfile('../images/fge_resized.tiff'):
    print('Resizing image...')
    subprocess.run([f'{venvName}/bin/python', 'resizeImage.py'], check=True)
