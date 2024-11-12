import subprocess
import os 

venvName='imageVenv'

if not os.path.isdir(venvName):
    subprocess.call(f'python3 -m venv {venvName}', shell=True)
    subprocess.call(f"bash -ic 'source {venvName}/bin/activate; pip3 install -r requirements.txt'", shell=True)

