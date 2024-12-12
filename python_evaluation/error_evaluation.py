import splinepy as sp
import numpy as np
import xml.etree.ElementTree as ET
from os import path
from splinepy.helpme.integrate import _get_integral_measure, _get_quadrature_information

DEGREE_ELEVATIONS = 1
WDIR = "--------------------------------------"
GEOMETRY_FILE = "-----------------------------"
PRESSURE_FILE = path.join(WDIR, "pressure_field_patches.xml")
VELOCITY_FILE = path.join(WDIR, "velocity_field_patches.xml")
PRESSURE_REC_FILE = path.join(WDIR, "pressure_field_rec_patches.xml")
VELOCITY_REC_FILE = path.join(WDIR, "velocity_field_rec_patches.xml")

common_show_options = {
    "cmap": "jet",
    "scalarbar": True,
    "lighting": "off",
    "control_points": False,
    "knots": False,
}

# Define error norm calculations
l1func = lambda orig, rec: np.abs(orig - rec)
l2func = lambda orig, rec: np.power(orig-rec, 2)
l1relfunc = lambda orig, rec: np.abs((orig-rec) / orig)
l2relfunc = lambda orig, rec: np.abs(np.power(orig-rec, 2) / orig)
norm_funcs = {
    "l1": l1func,
    "l2": l2func,
    "l1_rel": l1relfunc,
    "l2_rel": l2relfunc
}

def integrate(bezier, field):
    """Integrate field over Bezier patch"""
    meas = _get_integral_measure(bezier)
    quad_positions, quad_weights = _get_quadrature_information(
        bezier, orders=None
    )
    result = np.einsum(
        "i...,i,i->...",
        field.evaluate(quad_positions),
        meas(bezier, quad_positions),
        quad_weights,
        optimize=True
    )
    return result

def integrate_multipatch(multipatch, fields):
    """Integrate field over multipatch geometry"""
    integration_sum = 0.0
    for patch, field_patch in zip(multipatch.patches, fields):
        integration_sum += integrate(patch, field_patch)
    return integration_sum

def compute_integral_error(multipatch, fields_original, fields_recreated, norm="l2"):
    fields_error = []
    for patch, field_orig, field_recreated in zip(multipatch.patches, fields_original, fields_recreated):
        field_error = patch.copy()
        field_error.cps = norm_funcs[norm](field_orig, field_recreated)
        fields_error.append(field_error)
    error_integral = integrate_multipatch(multipatch, fields_error)
    if norm == "l2" or norm == "l2_rel":
        return np.sqrt(error_integral)
    else:
        return error_integral

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

def show_multipatch_field(
    mp,
    solution_vectors,
    data_name="solution",
    output_file=None,
    vmin=None,
    vmax=None
):
    """
    Show a field defined on a multipatch
    
    Parameters
    ------------
    mp: splinepy Multipatch
        Multipatch geometry object
    solution_vectors: list<np.ndarray>
        List of patch solution vectors
    data_name: str
        Name for field (e.g. velocity)
    """
    assert isinstance(solution_vectors, list), "Solution vectors have to be a list"
    assert len(mp.patches) == len(solution_vectors), "Mismatch between number of patches and patch solution vectors"
    spline_data_list = []
    for mp_patch, sv in zip(mp.patches, solution_vectors):
        data_patch = mp_patch.copy()
        data_patch.cps = sv
        spline_data_list.append(data_patch)
        mp_patch.spline_data[data_name] = spline_data_list[-1]
        mp_patch.show_options(data=data_name, **common_show_options)
    data_mp = sp.Multipatch(splines=spline_data_list)
    mp.spline_data[data_name] = data_mp
    mp.show_options(data=data_name, **common_show_options)
    if output_file is not None:
        assert isinstance(output_file, str), "Output file must be a string"
        assert vmin is not None and vmax is not None, "vmin and vmax must be given for scalarbar"
        sp.io.svg.export(
            output_file,
            mp.patches,
            scalarbar=True,
            vmin=vmin,
            vmax=vmax
        )
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
    pressure_rec_data = get_solution_vectors(file_path=PRESSURE_REC_FILE)
    velocity_rec_data = get_solution_vectors(file_path=VELOCITY_REC_FILE, two_dimensional=True)
    
    # Load the geometry
    microstructure, ms_vel = load_geometry(GEOMETRY_FILE, degree_elevations=DEGREE_ELEVATIONS)
    
    # Show pressure and velocity field
    show_multipatch_field(microstructure, pressure_data, data_name="pressure")
    show_multipatch_field(
        ms_vel,
        velocity_data,
        data_name="velocity",
        output_file="velocity_field.svg",
        vmin=0.0,           # Caution: vmin and vmax for scalarbar are hardcoded right now
        vmax=0.00119
    )
    
    # Compute errors
    norms = list(norm_funcs.keys())
    error_values = []
    print("-"*80)
    print("Pressure")
    for norm in norms:
        error_value = compute_integral_error(
            multipatch=microstructure,
            fields_original=pressure_data,
            fields_recreated=pressure_rec_data,
            norm=norm
        )
        error_values.append(error_value)
        print(norm, error_value)
        
    print("-"*80)
    print("Velocity")
    error_values = []
    for norm in norms:
        error_value = compute_integral_error(
            multipatch=ms_vel,
            fields_original=velocity_data,
            fields_recreated=velocity_rec_data,
            norm=norm
        )
        error_values.append(error_value)
        print(norm, error_value)