

import sys
import numpy as np
import yaml
from paraview.simple import *
from vtk.numpy_interface import dataset_adapter as dsa
from paraview.vtk.vtkFiltersGeneral import vtkGradientFilter
from paraview import servermanager
from paraview.simple import TrivialProducer

yml_file = sys.argv[1]
pvd_file = sys.argv[2]
output_file = sys.argv[3]

input_file = pvd_file

# ========== Load PVD File ==========
# Load data
reader = PVDReader(FileName=input_file)
reader.UpdatePipeline()
timesteps = reader.TimestepValues
print("Available timesteps:", timesteps)
# Choose last timestep
reader.UpdatePipeline(timesteps[-1])
# Use this to get the raw VTK data object in pvpython (✅ correct way)
output_data = reader.GetClientSideObject().GetOutputDataObject(0)
input0 = dsa.WrapDataObject(output_data)

# ========== Extract Data Arrays ==========
cminus = input0.PointData['cminus']
cplus = input0.PointData['cplus']
phi = input0.PointData['phi']

# ========== Constants ==========
mplus = 0.000264357526226176
mminus = 0.009868520605112863
m0 = 0.012523878125000

# ========== Compute Temporary Arrays ==========
ln_cplus = np.log(np.maximum(cplus * mplus, 1e-12))
input0.PointData.append(ln_cplus, 'ln_cplus_temp')

ln_sum_plus = -mplus / m0 * np.log(np.maximum(1 - mplus * cplus - mminus * cminus, 1e-12))
input0.PointData.append(ln_sum_plus, 'ln_sum_plus_temp')

ln_cminus = np.log(np.maximum(cminus * mminus, 1e-12))
input0.PointData.append(ln_cminus, 'ln_cminus_temp')

ln_sum_minus = -mminus / m0 * np.log(np.maximum(1 - mplus * cplus - mminus * cminus, 1e-12))
input0.PointData.append(ln_sum_minus, 'ln_sum_minus_temp')

# ========== Gradient Helper ==========
def compute_gradient(input_data, array_name, result_name):
    grad = vtkGradientFilter()
    grad.SetInputData(input_data.VTKObject)
    grad.SetInputArrayToProcess(0, 0, 0, 0, array_name)
    grad.SetResultArrayName(result_name)
    grad.Update()
    return dsa.WrapDataObject(grad.GetOutput()).PointData[result_name]

# ========== Compute Gradients ==========
grad_ln_cplus = compute_gradient(input0, 'ln_cplus_temp', 'grad_ln_cplus')
grad_ln_sum_plus = compute_gradient(input0, 'ln_sum_plus_temp', 'grad_ln_sum_plus')
grad_phi = compute_gradient(input0, 'phi', 'grad_phi')
grad_ln_cminus = compute_gradient(input0, 'ln_cminus_temp', 'grad_ln_cminus')
grad_ln_sum_minus = compute_gradient(input0, 'ln_sum_minus_temp', 'grad_ln_sum_minus')

# ========== Compute Fluxes ==========
input0.PointData.append(grad_ln_cplus * cplus, 'diff_flux_cplus')
input0.PointData.append(grad_ln_sum_plus * cplus, 'diff_size_cplus')
input0.PointData.append(grad_phi * cplus, 'migration_flux_cplus')
input0.PointData.append(cplus * (grad_ln_cplus + grad_ln_sum_plus + grad_phi), 'total_flux_cplus')

input0.PointData.append(grad_ln_cminus * cminus, 'diff_flux_cminus')
input0.PointData.append(grad_ln_sum_minus * cminus, 'diff_size_cminus')
input0.PointData.append(-grad_phi * cminus, 'migration_flux_cminus')
input0.PointData.append(cminus * (grad_ln_cminus + grad_ln_sum_minus - grad_phi), 'total_flux_cminus')

# ========== Save Output ==========
# Create a ParaView proxy from your modified VTK data
proxy = TrivialProducer()
proxy.GetClientSideObject().SetOutput(input0.VTKObject)

# Save the data using the proxy
SaveData(output_file, proxy=proxy)
print(" Flux computation complete. Output saved to flux_output.vtu.")
