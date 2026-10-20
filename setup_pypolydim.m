pyenv('Version', './.venv/bin/python','ExecutionMode','InProcess');
pyrunfile("pybind11_builtins.py");

pd = py.importlib.import_module('pypolydim');
export_vtk_utilities = py.importlib.import_module('pypolydim.export_vtk_utilities');
np = py.importlib.import_module('numpy');
polydim = pd.polydim;
gedim = pd.gedim;

% Rename namespace
MeshGenerator_Types_2D = pd.polydim.pde_tools.mesh.pde_mesh_utilities.MeshGenerator_Types_2D;
local_space_pcc_2_d = polydim.pde_tools.local_space_pcc_2_d;
local_space_pcc_3_d = polydim.pde_tools.local_space_pcc_3_d;
ExportVTKUtilities = export_vtk_utilities.ExportVTKUtilities;
DOFsData = polydim.pde_tools.do_fs.DOFsManager.DOFsData;
MeshGenerator_Types_3D = pd.polydim.pde_tools.mesh.pde_mesh_utilities.MeshGenerator_Types_3D;

% Overload function
m2py_float = @(x) py.numpy.array(x, pyargs('dtype', py.numpy.float64));
py2m_int =@(x) int64(py.numpy.array(x, pyargs('dtype', py.numpy.int64)));

export_solution_2 = @(file_name, mesh, cell0_ds_numeric, cell0_ds_exact, cell2_ds_error_l2, cell2_ds_error_h1) ...
    ExportVTKUtilities.export_solution_2(file_name, ...
        mesh, py.numpy.squeeze(py.numpy.array(cell0_ds_numeric, pyargs('dtype', py.numpy.float64))), ...
        py.numpy.squeeze(py.numpy.array(cell0_ds_exact, pyargs('dtype', py.numpy.float64))), ...
        py.numpy.squeeze(py.numpy.array(cell2_ds_error_l2, pyargs('dtype', py.numpy.float64))), ...
        py.numpy.squeeze(py.numpy.array(cell2_ds_error_h1, pyargs('dtype', py.numpy.float64))));

export_solution_3 = @(file_name, mesh, cell0_ds_numeric, cell0_ds_exact, cell3_ds_error_l2, cell3_ds_error_h1) ...
    ExportVTKUtilities.export_solution_3(file_name, ...
        mesh, py.numpy.squeeze(py.numpy.array(cell0_ds_numeric, pyargs('dtype', py.numpy.float64))), py.numpy.squeeze(py.numpy.array(cell0_ds_exact, pyargs('dtype', py.numpy.float64))), ...
        py.numpy.squeeze(py.numpy.array(cell3_ds_error_l2, pyargs('dtype', py.numpy.float64))), ...
        py.numpy.squeeze(py.numpy.array(cell3_ds_error_h1, pyargs('dtype', py.numpy.float64))));


% C++ Variable
MAX_UINT = 4294967295;

