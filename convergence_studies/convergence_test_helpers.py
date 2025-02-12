from os import chdir, getcwd
import numpy as np
import matplotlib.pyplot as plt
import subprocess
import re
from tqdm import tqdm
import scipy.sparse as scsp

cwd = getcwd()

def extract_number(keyword, text, number_pattern = "\d[-\dena.]*"):
    return re.findall(number_pattern, re.findall(f"{keyword}[\s]*{number_pattern}", text)[0])[0]

def find_error_values(text, relevant_pattern = "Errors:[\s]+Velocity: [-.\dena]*[\s]+Pressure: [-.\dena]*"):
    relevant_text = re.findall(relevant_pattern, text)[0]
    relevant_numbers = re.findall("[-.\dena]+", relevant_text.replace("Velocity", "Vlocity").replace("Pressure", "Prssur"))
    return (float(number) for number in relevant_numbers)

def extract_rhs(text, keyword="rhs"):
    relevant_pattern= f"{keyword}---[\s\dena.-]*---{keyword}"
    relevant_text = re.findall(relevant_pattern, text)[0][7:-6]
    return np.fromstring(relevant_text, sep="\n")

def plot_convergence_results(
    plottitle,
    velocity_orders,
    vel_errors_list,
    pressure_errors_list,
    ignore_pressure_error
):
    col_length = 6
    ncols = 1 if ignore_pressure_error else 2
    temporary_labels = ["vel. order = REPLACE"] + ["pres. order = REPLACE"] * (ncols-1)
    # xlabel = "No. of h-refinements"
    xlabel = "Element sizes"
    ylabels = [
        r"$|| \boldsymbol{v} - \boldsymbol{v}_{\text{ana}} ||_{L_2}$",
        r"$|| p - p_{\text{ana}} ||_{L_2}$"
    ]
    if ignore_pressure_error:
        ylabels.pop(1)
    
    fig,axes = plt.subplots(ncols=ncols, figsize=(ncols*col_length,4))
    if ncols == 1:
        axes = np.array([axes])
    for i,velocity_order in enumerate(velocity_orders):
        vel_errors = vel_errors_list[i,:]
        pressure_errors = pressure_errors_list[i,:]
        graph_values = [vel_errors]
        if not ignore_pressure_error:
            graph_values.append(pressure_errors)
            
        labels = temporary_labels.copy()
        element_sizes = np.power(2.0, -np.arange(len(vel_errors)))
        order_values = [velocity_order] + [velocity_order-1]*(ncols-1)
        labels = [text.replace("REPLACE", str(value)) for text,value in zip(labels, order_values)]
        
        for axis,graph_value,label in zip(axes, graph_values, labels):
            axis.loglog(element_sizes, graph_value, "^-", label=label)
        
    for axis,ylabel in zip(axes, ylabels):
        axis.set_xlabel(xlabel)
        axis.set_ylabel(ylabel)
        axis.set_aspect('equal')
        
    plt.legend()
    plt.title(plottitle)
    plt.axis("equal")
    plt.tight_layout()
    plt.show()
    
def run_convergence_study(executable_command_func, xml_file, n_refinements, p_refinements=[0,1,2,3]):
    """
    Runs a convergence study with a given executable of a (fluid) simulation. Plots the
    results in a loglog-plot
    
    Parameters
    ------------
    executable_command_list: callable -> list<str...>
        A function which returns of CLI commands to execute. This executable should have
        the following entries:
            -) href: number of h-refinements to perform
            -) pref: number of p-refinements to perform
    xml_file: str
        Filename of the simulation's xml-file
    n_refinements: int
        Number of h-refinements to perform
    p_refinements: list<int>
        List of which number of p-refinements to perform
    
    """
    velocity_orders = [2 + i for i in range(len(p_refinements))]   # using Taylor-Hood

    all_vel_errors_list = []
    all_p_errors_list = []

    # Go through each p-refinement
    for p_refinement in p_refinements:
        vel_errors_list = []
        pressure_errors_list = []
        # Go through each h-refinement
        for h_refinement  in tqdm(range(n_refinements)):
            # Execute simulation
            command_list = executable_command_func(href=h_refinement, pref=p_refinement)
            command_list.append(f"-f {xml_file}")
            program_output = subprocess.Popen(
                command_list,
                stdout=subprocess.PIPE,
                text=True
            ).communicate()[0]
            # Get simulation information
            velocity_error, pressure_error = find_error_values(program_output)
            vel_errors_list.append(velocity_error)
            pressure_errors_list.append(pressure_error)
        # Append all errors from current p-refinement to list of all errors
        all_vel_errors_list.append(vel_errors_list)
        all_p_errors_list.append(pressure_errors_list)
        
    # Convert error lists to numpy array
    all_vel_errors_list = np.array(all_vel_errors_list)
    all_p_errors_list = np.array(all_p_errors_list)

    plot_convergence_results(
        xml_file,
        velocity_orders,
        all_vel_errors_list,
        all_p_errors_list,
        ignore_pressure_error=False
    )