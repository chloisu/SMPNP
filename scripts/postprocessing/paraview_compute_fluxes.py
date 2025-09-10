import sys
import numpy as np
import yaml
import constants
from paraview.simple import *
from vtk.numpy_interface import dataset_adapter as dsa
from paraview.vtk.vtkFiltersGeneral import vtkGradientFilter
from paraview import servermanager
from paraview.simple import TrivialProducer

yml_file = sys.argv[1]
pvd_file = sys.argv[2]
output_file = sys.argv[3]

# ========== Get parameters from .yml file ==========
with open(yml_file, 'r') as f:
    mydict = yaml.safe_load(f)

# Example: extract values from the config
grid_type = mydict["grid_type"]
h = mydict["grid"][grid_type]["reservoir_height"]
r = mydict["grid"][grid_type]["pore_radius"]
l = mydict["grid"][grid_type]["pore_length"]
L_REF = mydict["non_dim"]["L_REF"]

# ========== Load PVD File ==========
# Load data
reader = PVDReader(FileName=pvd_file)
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
mplus = mydict["species_parameters"]["ion_sizes"][constants.CATION_EQ-2]**3*constants.MOL_PER_LITER_TO_PER_CUBIC_METER#0.000264357526226176
mminus = mydict["species_parameters"]["ion_sizes"][constants.ANION_EQ -2]**3*constants.MOL_PER_LITER_TO_PER_CUBIC_METER#0.009868520605112863
m0 = mydict["species_parameters"]["solvent_size"]**3*constants.MOL_PER_LITER_TO_PER_CUBIC_METER#0.012523878125000

# ========== Compute Temporary Arrays ==========
ln_cplus = np.log(np.maximum(cplus * mplus, 1e-12))
input0.PointData.append(ln_cplus, 'ln_cplus_temp')

ln_sum_plus = -mplus / m0 * np.log(np.maximum(1 - mplus * cplus - mminus * cminus, 1e-12))
input0.PointData.append(ln_sum_plus, 'ln_sum_plus_temp')

ln_cminus = np.log(np.maximum(cminus * mminus, 1e-12))
input0.PointData.append(ln_cminus, 'ln_cminus_temp')

ln_sum_minus = -mminus / m0 * np.log(np.maximum(1 - mplus * cplus - mminus * cminus, 1e-12))
input0.PointData.append(ln_sum_minus, 'ln_sum_minus_temp')

# ========== Chemical Potential ==========
mu_plus = np.log(cplus*mplus) - mplus/m0*np.log(np.maximum(1 - mplus * cplus - mminus * cminus, 1e-12)) + phi
input0.PointData.append(mu_plus, 'mu_plus')

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

# ========== Compute Dimensional Quantities ==========
# extract reference quantities
C_REF = float(mydict["non_dim"]["C_REF"])
PHI_REF = float(mydict["non_dim"]["PHI_REF"])
L_REF = float(mydict["non_dim"]["L_REF"])
D_REF = float(mydict["non_dim"]["D_REF"])

# extract other inputs
D_plus = float(mydict["species_parameters"]["diffusivities"][constants.CATION_EQ-2])
D_minus = float(mydict["species_parameters"]["ion_sizes"][constants.ANION_EQ -2])
# primary variables
cminus_1_per_cubic_meter = cminus * C_REF
cplus_1_per_cubic_meter= cplus * C_REF
phi_volt = phi * PHI_REF
KB_T = constants.K_B*constants.T
mu_plus_diff_joule = KB_T*ln_cplus
mu_plus_size_joule = KB_T*ln_sum_plus
mu_plus_potential_joule = constants.E_CHARGE*PHI_REF*phi

mu_minus_diff_joule = KB_T*ln_cminus
mu_minus_size_joule = KB_T*ln_sum_minus
mu_minus_potential_joule = -constants.E_CHARGE*PHI_REF*phi
## fluxes
# diffusion
diff_flux_cplus_1_per_meter_squared_per_second = D_REF*C_REF/L_REF*(D_plus/D_REF)*cplus*grad_ln_cplus
diff_flux_cminus_1_per_meter_squared_per_second = D_REF*C_REF/L_REF*(D_minus/D_REF)*cminus*grad_ln_cminus
# migration
migration_flux_cplus_1_per_meter_squared_per_second = D_REF*C_REF/L_REF*(D_plus/D_REF)*cplus*grad_phi
migration_flux_cminus_1_per_meter_squared_per_second = D_REF*C_REF/L_REF*(D_minus/D_REF)*cminus*(-1.0)*grad_phi
# size
size_flux_cplus_1_per_meter_squared_per_second =  D_REF*C_REF/L_REF*(D_plus/D_REF)*cplus*grad_ln_sum_plus
size_flux_cminus_1_per_meter_squared_per_second = D_REF*C_REF/L_REF*(D_minus/D_REF)*cminus*grad_ln_sum_minus
# total
tota_flux_cplus_1_per_meter_squared_per_second = diff_flux_cplus_1_per_meter_squared_per_second + migration_flux_cplus_1_per_meter_squared_per_second + size_flux_cplus_1_per_meter_squared_per_second
tota_flux_cminus_1_per_meter_squared_per_second = diff_flux_cminus_1_per_meter_squared_per_second + migration_flux_cminus_1_per_meter_squared_per_second + size_flux_cminus_1_per_meter_squared_per_second
## output
input0.PointData.append(cplus_1_per_cubic_meter, 'cplus_1_per_cubic_meter')
input0.PointData.append(cminus_1_per_cubic_meter, 'cminus_1_per_cubic_meter')
input0.PointData.append(phi_volt, 'phi_volt')

input0.PointData.append(mu_plus_diff_joule, 'mu_plus_diff_joule')
input0.PointData.append(mu_plus_size_joule, 'mu_plus_size_joule')
input0.PointData.append(mu_plus_potential_joule, 'mu_plus_potential_joule')

input0.PointData.append(mu_minus_diff_joule, 'mu_minus_diff_joule')
input0.PointData.append(mu_minus_size_joule, 'mu_minus_size_joule')
input0.PointData.append(mu_minus_potential_joule, 'mu_minus_potential_joule')

input0.PointData.append(diff_flux_cplus_1_per_meter_squared_per_second, 'diff_flux_cplus_1_per_meter_squared_per_second')
input0.PointData.append(diff_flux_cminus_1_per_meter_squared_per_second, 'diff_flux_cminus_1_per_meter_squared_per_second')

input0.PointData.append(migration_flux_cplus_1_per_meter_squared_per_second, 'migration_flux_cplus_1_per_meter_squared_per_second')
input0.PointData.append(migration_flux_cminus_1_per_meter_squared_per_second, 'migration_flux_cminus_1_per_meter_squared_per_second')

input0.PointData.append(size_flux_cplus_1_per_meter_squared_per_second, 'size_flux_cplus_1_per_meter_squared_per_second')
input0.PointData.append(size_flux_cminus_1_per_meter_squared_per_second, 'size_flux_cminus_1_per_meter_squared_per_second')

input0.PointData.append(tota_flux_cplus_1_per_meter_squared_per_second, 'tota_flux_cplus_1_per_meter_squared_per_second')
input0.PointData.append(tota_flux_cminus_1_per_meter_squared_per_second, 'tota_flux_cminus_1_per_meter_squared_per_second')

# ========== Save Output ==========
# Create a ParaView proxy from your modified VTK data
proxy = TrivialProducer()
proxy.GetClientSideObject().SetOutput(input0.VTKObject)

# Save the data using the proxy
SaveData(output_file, proxy=proxy)
print(" Flux computation complete.")
