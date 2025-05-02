import numpy as np
from vtk.numpy_interface import dataset_adapter as dsa
from vtk.numpy_interface import algorithms as algs
from paraview.vtk.vtkFiltersGeneral import vtkGradientFilter

# Wrap the input
input0 = dsa.WrapDataObject(self.GetInput())

# Get raw data arrays
cminus = input0.PointData['cminus']
cplus = input0.PointData['cplus']
phi = input0.PointData['phi']

# define the required mvalues
mplus = 0.000264357526226176
mminus = 0.009868520605112863
m0 = 0.012523878125000
# Compute ln(cplus), safely
ln_cplus = np.log(np.maximum(cplus*mplus, 1e-12))  # avoid log(0)
input0.PointData.append(ln_cplus, 'ln_cplus_temp')

# Compute ln(1-sum) safely
ln_sum_plus = -mplus/m0*np.log(np.maximum((1-mplus*cplus-mminus*cminus),1e-12))
input0.PointData.append(ln_sum_plus, 'ln_sum_plus_temp')

# Compute ln(cminus), safely
ln_cminus = np.log(np.maximum(cminus*mminus, 1e-12))  # avoid log(0)
input0.PointData.append(ln_cminus, 'ln_cminus_temp')

# Compute ln(1-sum) safely
ln_sum_minus = -mminus/m0*np.log(np.maximum((1-mplus*cplus-mminus*cminus),1e-12))
input0.PointData.append(ln_sum_minus, 'ln_sum_minus_temp')


# Function to compute gradient via vtkGradientFilter
def compute_gradient(input_data, array_name, result_name):
    grad = vtkGradientFilter()
    grad.SetInputData(input_data)
    grad.SetInputArrayToProcess(0, 0, 0, 0, array_name)
    grad.SetResultArrayName(result_name)
    grad.Update()
    return dsa.WrapDataObject(grad.GetOutput()).PointData[result_name]

# Compute gradients
grad_ln_cplus = compute_gradient(self.GetInput(), 'ln_cplus_temp', 'grad_ln_cplus')
grad_ln_sum_plus = compute_gradient(self.GetInput(), 'ln_sum_plus_temp', 'grad_ln_sum_minus')
grad_phi = compute_gradient(self.GetInput(), 'phi', 'grad_phi')

# Compute gradients
grad_ln_cminus = compute_gradient(self.GetInput(), 'ln_cminus_temp', 'grad_ln_cminus')
grad_ln_sum_minus = compute_gradient(self.GetInput(), 'ln_sum_minus_temp', 'grad_ln_sum_plus')

# Compute fluxes
output = dsa.WrapDataObject(self.GetOutput())
output.PointData.append(grad_ln_cplus*cplus, 'diff_flux_cplus')
output.PointData.append(grad_ln_sum_plus*cplus, 'diff_size_cplus')
output.PointData.append(grad_phi*cplus, 'migration_flux_cplus')
output.PointData.append(cplus*(grad_ln_cplus + grad_ln_sum_plus+grad_phi), 'total_flux_cplus')

output.PointData.append(grad_ln_cminus*cminus, 'diff_flux_cminus')
output.PointData.append(grad_ln_sum_minus*cminus, 'diff_size_cminus')
output.PointData.append(-grad_phi*cminus, 'migration_flux_cminus')
output.PointData.append(cminus*(grad_ln_cminus + grad_ln_sum_minus-grad_phi), 'total_flux_cminus')
