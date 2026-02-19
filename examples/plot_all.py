import glob
import os
import subprocess

cwd = os.path.dirname(os.path.realpath(__file__))
dirs = glob.glob(os.path.join(cwd,'*alert/'))
dirs.sort()

for d in dirs:
    name_dir = d.split("/")[-2]
    print(name_dir)
    plot_script = glob.glob(os.path.join(d, "plot_*"))[0]
    subprocess.run(
        f'python {plot_script}',
        shell=True,
        executable="/bin/bash"
    )