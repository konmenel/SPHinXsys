import os
import sys
import re
from pathlib import Path
from paraview.simple import *


PATH = Path(os.getenv('HOME')) / "SPHinXsys-build"/"cases_user"/"2d_dambreak_pinned_boulder"
if sys.platform.startswith('linux'):
    PATH = str(PATH) + "/bin/output/"
else:
    PATH = str(PATH) + "\\src\\output\\"


def open_files():
    wall_re = re.compile(r'SPHBody_Wall_[0-9]+.vtu')
    boulder_re = re.compile(r'SPHBody_Boulder_[0-9]+.vtu')
    water_re = re.compile(r'SPHBody_WaterBody_[0-9]+.vtu')
    wall_files = []
    boulder_files = []
    water_files = []
    
    files = os.listdir(PATH)
    for name in files:
        if boulder_re.fullmatch(name):
            boulder_files.append(PATH + name)
        elif water_re.fullmatch(name):
            water_files.append(PATH + name)
        elif wall_re.fullmatch(name):
            wall_files.append(PATH + name)
    
    num_re = re.compile(r'_([0-9]+).vtu$')
    for list_ in [water_files, boulder_files]:
        list_.sort(key=lambda x: int(num_re.findall(x)[0]))

    return (OpenDataFile(wall_files), OpenDataFile(boulder_files),
            OpenDataFile(water_files), max((len(water_files)), len(boulder_files)))


ResetSession()
# for x in GetSources().values():
#     Delete(x[0])

view = GetActiveView()
wall, boulder, water, number_of_files = open_files()
scene = GetAnimationScene()
scene.PlayMode = 'Snap To TimeSteps'

RenameSource('Wall', wall)
RenameSource('Boulder', boulder)
RenameSource('WaterBody', water)

wall_disp = GetDisplayProperties(wall, view)
boulder_disp = GetDisplayProperties(boulder, view)
water_disp = GetDisplayProperties(water, view)

for disp in (wall_disp, boulder_disp, water_disp):
    disp.SetRepresentationType('Point Gaussian')
    disp.GaussianRadius = 0.0017

boulder_disp.DiffuseColor = [0.67, 0.33, 0.0]

for source in [wall, boulder]:
    Show(source, view)

display_w = Show(water)

ColorBy(display_w, ('POINTS', 'Velocity'), separate = True)
water_color_func = GetColorTransferFunction('Velocity', display_w, separate=True)
water_color_func.ApplyPreset('Jet', True)

scene.GoToLast()
scene.GoToPrevious()
water_color_func.RescaleTransferFunctionToDataRange()
scene.GoToFirst()

display_w.SetScalarBarVisibility(view, True)

Render()
