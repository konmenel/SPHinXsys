import os
import sys
import re
from pathlib import Path
from paraview.simple import *


PARTICLE_RADIUS = 0.0013 / 2. / 2.

PATH = Path(os.getenv('HOME')) / "SPHinXsys-build"/"tests"/"user_examples"/"2d_boulder_wall"
if sys.platform.startswith('linux'):
    PATH = str(PATH) + "/bin/output/"
else:
    PATH = str(PATH) + "\\src\\output\\"


def open_files():
    wall_re = re.compile(r'SPHBody_Wall_[0-9.]+.vtp$')
    boulder_re = re.compile(r'SPHBody_Boulder_[0-9.]+.vtp$')
    wall_files = []
    boulder_files = []
    
    files = os.listdir(PATH)
    for name in files:
        if boulder_re.fullmatch(name):
            boulder_files.append(PATH + name)
        elif wall_re.fullmatch(name):
            wall_files.append(PATH + name)
    
    num_re = re.compile(r'_([0-9]+).vtp$')
    boulder_files.sort(key=lambda x: int(num_re.findall(x)[0]))

    return (OpenDataFile(wall_files), OpenDataFile(boulder_files),
            len(boulder_files))


ResetSession()
# for x in GetSources().values():
#     Delete(x[0])

view = GetActiveView()
wall, boulder, number_of_files = open_files()
scene = GetAnimationScene()
scene.PlayMode = 'Snap To TimeSteps'

RenameSource('Wall', wall)
RenameSource('Boulder', boulder)

wall_disp = GetDisplayProperties(wall, view)
boulder_disp = GetDisplayProperties(boulder, view)

for disp in (wall_disp, boulder_disp):
    disp.SetRepresentationType('Point Gaussian')
    disp.GaussianRadius = PARTICLE_RADIUS

boulder_disp.DiffuseColor = [0.67, 0.33, 0.0]

for source in [wall, boulder]:
    Show(source, view)

Render()