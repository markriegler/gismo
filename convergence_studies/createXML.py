import splinepy as sp
from export_helpers import export, AdditionalBlocks
import json
from os import path

CASES_FILE = "stokes_studies.json"

ASSEMBLY_OPTIONS_ID = 10
VEL_BC_ID = 1
P_BC_ID = 2
BODY_FORCE_ID = 100
VEL_ANALYTICAL_ID = 101
P_ANALYTICAL_ID = 102
PRESSURE_ID = 1
VELOCITY_ID = 0

def get_case_dict(casename):
    with open(CASES_FILE, "r") as f:
        studies_dict = json.load(f)
    return studies_dict[casename]

def bottom_identifier(points):
    return points[:,1] < 1e-8

def create_xml(casename, outdir):
    case_dict = get_case_dict(casename)
    
    outname = path.join(outdir, f"{case_dict['outname']}.xml")
        
    # Process the geometry type
    geometry_type = case_dict["geometry"]
    if geometry_type == "box2D":
        geometry = sp.helpme.create.box(1,1)
    elif geometry_type == "box3D":
        geometry = sp.helpme.create.box(1,1,1)
    else:
        raise ValueError(f"Cannot process geometry type {geometry_type}")
        
    additional_blocks = AdditionalBlocks()
    
    additional_blocks.add_assembly_options(block_id=ASSEMBLY_OPTIONS_ID, comment=" Assembly options ", dirichlet_values=102)
    
    vel_analytical = tuple(function_string for function_string in case_dict["vel_analytical"].values())
    pres_analytical = case_dict["p_analytical"]
    
    if "vel_bcs" in case_dict.keys():
        vel_bcs = [tuple(function_string for function_string in case_dict["vel_bcs"].values())]
    else:
        vel_bcs = [vel_analytical]
    
    additional_blocks.add_boundary_conditions(
        block_id=VEL_BC_ID,
        dim=geometry.dim,
        function_list=vel_bcs,
        bc_list=[("BID1", "Dirichlet", 0), ("BID2", "Dirichlet", 0)],
        unknown_id=VELOCITY_ID,
        comment=" Velocity boundary conditions "
    )
    
    # For pressure boundary conditions just apply corner value
    if "p_se_corner" in case_dict.keys():
        p_se_corner_value = case_dict["p_se_corner"]
    else:
        p_se_corner_value = 0.0
    additional_blocks.add_boundary_conditions(
        block_id=P_BC_ID,
        dim=geometry.dim,
        cv_list=[
            (PRESSURE_ID, 0, 1, p_se_corner_value)
        ],
        unknown_id=PRESSURE_ID,
        comment=" Pressure boundary conditions "
    )
    
    additional_blocks.add_function(
        dim=geometry.dim,
        block_id=BODY_FORCE_ID,
        function_string=tuple(function_string for function_string in case_dict["bodyforce"].values()),
        comment=" Body forces "
    )
    
    additional_blocks.add_function(
        dim=geometry.dim,
        block_id=VEL_ANALYTICAL_ID,
        function_string=vel_analytical,
        comment=" Analytical velocity "
    )
    
    additional_blocks.add_function(
        dim=geometry.dim,
        block_id=P_ANALYTICAL_ID,
        function_string=pres_analytical,
        comment=" Analytical pressure "
    )
    
    geometry = sp.Multipatch(splines=[geometry])
    if geometry_type == "box2D":
        geometry.boundary_from_function(function=bottom_identifier, boundary_id=2)
    
    export(
        fname=outname,
        multipatch=geometry,
        additional_blocks=additional_blocks.to_list()
    )
    
    print(f"Exported xml-file to {outname}")
    
if __name__ == "__main__":
    for casename in ["moeller", "buffa", "vortexflux"]:
        create_xml(casename, ".")