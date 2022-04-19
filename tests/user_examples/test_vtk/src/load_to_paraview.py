import os
import sys
import re
from pathlib import Path
# from paraview.simple import *

parent_dirname = Path(__file__).parent.parent.name
PATH = Path(os.getenv('HOME')) / "SPHinXsys-build"/"tests"/"user_examples"/parent_dirname
if sys.platform.startswith('linux'):
    PATH = str(PATH) + "/bin/output/"
else:
    PATH = str(PATH) + "\\src\\output\\"


def open_files():
    wall_re = re.compile(r'SPHBody_Wall_[0-9]+.vtk$')
    water_re = re.compile(r'SPHBody_WaterBody_[0-9]+.vtk$')
    wall_files = []
    water_files = []
    
    files = os.listdir(PATH)
    for name in files:
        if water_re.fullmatch(name):
            water_files.append(PATH + name)
        elif wall_re.fullmatch(name):
            wall_files.append(PATH + name)
    
    num_re = re.compile(r'_([0-9]+).vtk$')
    water_files.sort(key=lambda x: int(num_re.findall(x)[0]))

    return (OpenDataFile(wall_files), OpenDataFile(water_files), 
            len(water_files))


ResetSession()
# for x in GetSources().values():
#     Delete(x[0])

view = GetActiveView()
wall, water, number_of_files = open_files()
scene = GetAnimationScene()
scene.PlayMode = 'Snap To TimeSteps'

RenameSource('Wall', wall)
RenameSource('WaterBody', water)

wall_disp = GetDisplayProperties(wall, view)
water_disp = GetDisplayProperties(water, view)

for disp in (wall_disp, water_disp):
    disp.SetRepresentationType('Point Gaussian')
    disp.GaussianRadius = 0.005

Show(wall, view)

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
