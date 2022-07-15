#!/usr/bin/env pvpython
"""This is a script to load a bunch of files into ParaView for visualization.
By default, it will load all the files in the current directory and all its
subdirectories, but optionally you can specify a directory.

The script will open a ParaView session and load all the vtk files found in the
root directory and all its subdirectories. It could be run from the terminal or
from the ParaView GUI `Run Script` button in the python shell (the latter is 
recommended). If run from the terminal, it should be run using `pvpython`.

This script is a template, you can customize the visualization by changing the
`add_visualization_option` function.
"""
import os
import re
from typing import Dict, List, Optional

from paraview.simple import *


def get_vtk_files(path: Optional[str]=None) -> Dict[str, List[str]]:
    """Get the vtk files in the current directory and all its subdirectories.

    Parameters
    ----------
    path : str
        The path of the directory to search.

    Returns
    -------
    Dict[str, List[str]]
        A dictionary with the key the name of the temporal series and the value a list of
        the vtk files found with the same name. 
        Example: {'data/cube': ['data/cube1.vtk', 'data/cube2.vtk']}
    """
    if path is None:
        path, _ = os.path.split(os.path.abspath(__file__))
    
    print(path)

    vtk_files = {}

    # The regular expression to match the vtk files.
    vtk_re = re.compile(r'(?P<filename>[A-Za-z\-_.]+)(?P<number>[\d]*).vt[kpu]$')

    # Walk through the directory and find the vtk files.    
    for root, _, files in os.walk(path):
        for file in files:
            match = vtk_re.match(file)
            if match:
                filename = match.group('filename')

                fullname = os.path.join(root, filename)

                if fullname not in vtk_files:
                    vtk_files[fullname] = []
                
                vtk_files[fullname].append(os.path.join(root, file))


    # Sort the files by number if the number exists.
    for files in vtk_files.values():
        if vtk_re.findall(files[0])[0][1]:
            files.sort(key=lambda x: int(vtk_re.findall(x)[0][1]))
    
    return vtk_files


def open_files(files_dict: Dict[str, List[str]]) -> dict:
    """Function to open the files in ParaView.

    Parameters
    ----------
    files_dict : Dict[str, List[str]]
        A dictionary with the key the name of the temporal series and the value a list of
        the vtk files.
        Example: {'data/cube': ['data/cube1.vtk', 'data/cube2.vtk']}

    Returns
    -------
    dict
        A dictionary with the key the name of the temporal series and the value the object
        created by ParaView's reader.
        Example for the input example above: {'cube': <instance of Paraview vtkReader>}
    """
    pvobjects_dict = {}
    for key, value in files_dict.items():
        _, name = os.path.split(key)
        reader = OpenDataFile(value)

        RenameSource(name, reader)

        pvobjects_dict[name] = reader

    return pvobjects_dict


def add_visualization_option(pvobjects_dict: dict) -> None:
    """Generic function that add visualization option to paraview

    Parameters
    ----------
    pvobjects_dict : dict
        A dictionary with the key the name of the temporal series and the value the object
        created by ParaView's reader.
        Example: {'cube': <instance of Paraview vtkReader>}
    """
    
    wall_re = re.compile(r'wall', flags=re.IGNORECASE)
    water_re = re.compile(r'water', flags=re.IGNORECASE)
    boulder_re = re.compile(r'boulder', flags=re.IGNORECASE)

    view = GetActiveViewOrCreate('RenderView')
    scene = GetAnimationScene()
    scene.PlayMode = 'Snap To TimeSteps'

    for name, pvobj in pvobjects_dict.items():
        rep = GetDisplayProperties(pvobj, view)
        # rep.SetRepresentationType('Point Gaussian')
        # rep.GaussianRadius = 0.0045

        if wall_re.search(name):
            ColorBy(rep, ('POINTS', ''))
            rep.Opacity = 0.1
            
        elif water_re.search(name):
            ColorBy(rep, ('POINTS', 'Velocity'), separate=True)
            water_color_func = GetColorTransferFunction('Velocity', rep, separate=True)
            water_color_func.ApplyPreset('Jet', True)
            scene.GoToLast()
            scene.GoToPrevious()
            water_color_func.RescaleTransferFunctionToDataRange()
            scene.GoToFirst()
            rep.SetScalarBarVisibility(view, True)

        elif boulder_re.search(name):
            ColorBy(rep, ('POINTS', ''))
            rep.DiffuseColor = (0.4, 0.219, 0.043)


def main() -> int:
    """Main function.

    Returns
    -------
    int
        Exit code.
    """
    vtk_files = get_vtk_files()
    pvobjects_dict = open_files(vtk_files)
    
    view = GetActiveViewOrCreate('RenderView')

    for source in pvobjects_dict.values():
        Show(source, view)

    add_visualization_option(pvobjects_dict)

    Render()

    return 0



if __name__ == "__main__":
    raise SystemExit(main())

elif __name__ == "__vtkconsole__":
    main()
