import splinepy as sp
import numpy as np
import xml.etree.ElementTree as ET
from os import path

DEGREE_ELEVATIONS = 1
WDIR = "where/is/it"
GEOMETRY_FILE = "where/is/it"
PRESSURE_FILE = path.join(WDIR, "pressure_field.xml")
VELOCITY_FILE = path.join(WDIR, "velocity_field.xml")

common_show_options = {
    "cmap": "jet",
    "scalarbar": True,
    "lighting": "off",
    "control_points": False,
    "knots": False,
}

# Load XML file to numpy array
def get_solution_vectors(file_path, two_dimensional=False):
    # Parse the XML file
    tree = ET.parse(file_path)
    root = tree.getroot()
    
    solution_vectors = []
    
    for child in root:
        solution_vector = np.fromstring(child.text.strip(), sep="\n")
        solution_vector = np.nan_to_num(solution_vector, copy=False, nan=0.0, posinf=0.0, neginf=0.0)
        if two_dimensional:
            solution_vector = solution_vector.reshape(-1, 2)
        else:
            solution_vector = solution_vector.reshape(-1, 1)
        solution_vectors.append(solution_vector)

    return solution_vectors

def show_multipatch_field(mp, solution_vectors, data_name="solution"):
    """
    Show a field defined on a multipatch
    
    Parameters
    ------------
    mp: splinepy Multipatch
        Multipatch geometry object
    solution_vectors: list<np.ndarray>
        List of patch solution vectors
    """
    assert isinstance(solution_vectors, list), "Solution vectors have to be a list"
    assert len(mp.patches) == len(solution_vectors), "Mismatch between number of patches and patch solution vectors"
    spline_data_list = []
    for mp_patch, sv in zip(mp.patches, solution_vectors):
        data_patch = mp_patch.copy()
        data_patch.cps = sv
        spline_data_list.append(data_patch)
        mp_patch.spline_data[data_name] = spline_data_list[-1]
    data_mp = sp.Multipatch(splines=spline_data_list)
    mp.spline_data[data_name] = data_mp
    mp.show_options(data=data_name, **common_show_options)
    mp.show()
    
def load_geometry(filename, degree_elevations=0):
    """
    Load geometry and for velocity perform one degree elevation (Taylor-Hood elements)
    
    filename: str
        Filename of xml-file
    """
    microstructure = sp.io.gismo.load(filename)[0]
    if degree_elevations > 0:
        [patch.elevate_degrees([0,1]*degree_elevations) for patch in microstructure.patches]
    patches_p_elevated = []
    for patch in microstructure.patches:
        patch_elevated = patch.copy()
        patch_elevated.elevate_degrees([0,1])
        patches_p_elevated.append(patch_elevated)
    microstructure_vel = sp.Multipatch(splines=patches_p_elevated)
    microstructure_vel.determine_interfaces()
    
    return microstructure, microstructure_vel


if __name__ == "__main__":
    # Get solution vector for pressure and velocity
    pressure_data = get_solution_vectors(file_path=PRESSURE_FILE)
    velocity_data = get_solution_vectors(file_path=VELOCITY_FILE, two_dimensional=True)
    
    # Load the geometry
    microstructure, ms_vel = load_geometry(GEOMETRY_FILE, degree_elevations=DEGREE_ELEVATIONS)
    
    # Show pressure and velocity field
    show_multipatch_field(microstructure, pressure_data, data_name="pressure")
    show_multipatch_field(ms_vel, velocity_data, data_name="velocity")