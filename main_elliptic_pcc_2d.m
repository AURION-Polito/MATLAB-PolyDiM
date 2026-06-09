%%
clc
clear
close all

addpath('./Elliptic_PCC_2D/')

setup_pypolydim

%profile on

%% Arguments

[~, ~, ~] = mkdir('./Export');
export_path = './Export/Elliptic_PCC_2D';
[~, ~, ~] = mkdir(export_path);

mesh_type = pd.polydim.pde_tools.mesh.pde_mesh_utilities.MeshGenerator_Types_2D.triangular;
method_type = polydim.pde_tools.local_space_pcc_2_d.MethodTypes.vem_pcc;
method_order = py.int(2);
import_path = './'; % Mesh import path

geometry_utilities_config = gedim.GeometryUtilitiesConfig();
geometry_utilities_config.tolerance1_d = 1.0e-12;
geometry_utilities_config.tolerance2_d = 1.0e-14;
geometry_utilities = gedim.GeometryUtilities(geometry_utilities_config);

mesh_utilities = gedim.MeshUtilities();
vtk_utilities = export_vtk_utilities.ExportVTKUtilities();

%% PDE Data
test_polynomial

%% Convergence test

max_relative_areas = [0.1, 0.05, 0.01, 0.005];
errors = [];


for ref=1:length(max_relative_areas)
    
    %% Create Mesh

    max_relative_area = py.float(max_relative_areas(ref));
    [mesh, mesh_data, mesh_geometric_data] = create_mesh(gedim, polydim, geometry_utilities, mesh_utilities, mesh_type, max_relative_area, import_path, pde_domain);

    export_mesh_path = export_path + "/Mesh" + int2str(ref);
    [~, ~, ~] = mkdir(export_mesh_path);
    vtk_utilities.export_mesh(export_mesh_path, mesh)


    %% Discrete Local Space

    reference_element_data = polydim.pde_tools.local_space_pcc_2_d.create_reference_element(method_type, method_order);
    mesh_connectivity_data = polydim.pde_tools.mesh.MeshMatricesDAO_mesh_connectivity_data(mesh);
    dof_manager = polydim.pde_tools.do_fs.DOFsManager();
    mesh_do_fs_info = polydim.pde_tools.local_space_pcc_2_d.set_mesh_do_fs_info(reference_element_data, mesh, boundary_info);
    do_fs_data = dof_manager.create_do_fs_2_d(mesh_do_fs_info, mesh_connectivity_data);
    do_fs_data_indices = dof_manager.compute_cells_do_fs_indices(do_fs_data, py.int(2));

    %% Assemble

    num_cell2D = int64(mesh.cell2_d_total_number());
    num_dofs = int64(do_fs_data.number_do_fs);
    num_strongs = int64(do_fs_data.number_strongs);
    cell_dimension = py.int(2);

    global_rhs = zeros(num_dofs, 1);
    solution_dirichlet = zeros(num_strongs, 1);

    global_lhs_values = cell(num_cell2D, 1);
    global_lhs_rows = cell(num_cell2D, 1);
    global_lhs_cols = cell(num_cell2D, 1);

    dirichlet_lhs_values = cell(num_cell2D, 1);
    dirichlet_lhs_rows = cell(num_cell2D, 1);
    dirichlet_lhs_cols = cell(num_cell2D, 1);

    equation = polydim.pde_tools.equations.EllipticEquation();


    for c = 1:num_cell2D

        cell_index = py.int(c - 1);

        local_space_data = polydim.pde_tools.local_space_pcc_2_d.create_local_space(geometry_utilities_config.tolerance1_d, geometry_utilities_config.tolerance2_d, mesh_geometric_data, cell_index, reference_element_data);

        basis_functions_values = polydim.pde_tools.local_space_pcc_2_d.basis_functions_values(reference_element_data, local_space_data);
        basis_functions_derivative_values = polydim.pde_tools.local_space_pcc_2_d.basis_functions_derivative_values(reference_element_data, local_space_data);

        cell2_d_internal_quadrature = polydim.pde_tools.local_space_pcc_2_d.internal_quadrature(reference_element_data, local_space_data);
        weights = cell2_d_internal_quadrature.weights;
        points = double(cell2_d_internal_quadrature.points);

        diffusion_term_values = diffusion_term(points);
        k_max = max(abs(diffusion_term_values));

        diffusion_term_values = py.numpy.array(diffusion_term_values, pyargs('dtype', py.numpy.float64));
        source_term_values = py.numpy.array(source_term(points), pyargs('dtype', py.numpy.float64));
        
        % Compute local stiffness matrix
        local_a = equation.compute_cell_diffusion_matrix(diffusion_term_values, basis_functions_derivative_values, weights);
        local_a = local_a + k_max * polydim.pde_tools.local_space_pcc_2_d.stabilization_matrix(reference_element_data, local_space_data);
        local_a = double(local_a);
        
        % Compute local rhs
        local_rhs = double(equation.compute_cell_forcing_term(source_term_values, basis_functions_values, weights))';

        global_do_fs = do_fs_data.cells_global_do_fs{cell_dimension}{cell_index};

        % Test slicing indices
        global_dof_index = int64(py.numpy.array(do_fs_data_indices.cells_do_fs_global_index{cell_index}, pyargs('dtype', py.numpy.int64)))' + 1;
        local_dof_index = int64(py.numpy.array(do_fs_data_indices.cells_do_fs_local_index{cell_index}, pyargs('dtype', py.numpy.int64)))' + 1;

        global_rhs(global_dof_index) = global_rhs(global_dof_index) + local_rhs(local_dof_index);

        % Trial slicing indices
        global_strong_index = int64(py.numpy.array(do_fs_data_indices.cells_strongs_global_index{cell_index}, pyargs('dtype', py.numpy.int64)))' + 1;
        local_strong_index = int64(py.numpy.array(do_fs_data_indices.cells_strongs_local_index{cell_index}, pyargs('dtype', py.numpy.int64)))' + 1;

        local_lhs_dofs = local_a(local_dof_index, local_dof_index);
        global_lhs_values{c} = local_lhs_dofs(:);
        global_lhs_rows{c} = repmat(global_dof_index, length(global_dof_index), 1);
        global_lhs_cols{c} = repelem(global_dof_index, length(global_dof_index), 1);

        local_lhs_strongs = local_a(local_dof_index, local_strong_index);
        dirichlet_lhs_values{c} = local_lhs_strongs(:);
        dirichlet_lhs_rows{c} = repmat(global_dof_index, length(global_strong_index), 1);
        dirichlet_lhs_cols{c} = repelem(global_strong_index, length(global_dof_index), 1);
        
        % Compute strong-dirichlet term
        solution_dirichlet = compute_strong_term(polydim, cell_index, mesh, mesh_do_fs_info, do_fs_data, reference_element_data, local_space_data, strong_boundary_condition, solution_dirichlet);
        
        % Add neumann boundary condition
        global_rhs = compute_weak_term(polydim, gedim, cell_index, mesh, mesh_geometric_data, mesh_do_fs_info, do_fs_data, reference_element_data, local_space_data, weak_boundary_condition, global_rhs);
    end

    global_matrix_a = sparse(vertcat(global_lhs_rows{:}), vertcat(global_lhs_cols{:}), vertcat(global_lhs_values{:}), num_dofs, num_dofs);
    dirichlet_matrix_a = sparse(vertcat(dirichlet_lhs_rows{:}), vertcat(dirichlet_lhs_cols{:}), vertcat(dirichlet_lhs_values{:}), num_dofs, num_strongs);

    global_rhs = global_rhs - dirichlet_matrix_a * solution_dirichlet;
    solution = global_matrix_a \ global_rhs;

    %% Compute Errors


    residual_norm = 0.0;

    if do_fs_data.number_do_fs > 0
        residual = global_matrix_a * solution - global_rhs;
        residual_norm = norm(residual);
    end

    numCell0D = int64(mesh.cell0_d_total_number());
    cell0_ds_coordinates = double(mesh.cell0_ds_coordinates());

    cell2_ds_error_l2 = zeros(num_cell2D, 1);
    cell2_ds_norm_l2 = zeros(num_cell2D, 1);
    cell2_ds_error_h1 = zeros(num_cell2D, 1);
    cell2_ds_norm_h1 = zeros(num_cell2D, 1);
    mesh_size = 0.0;

    for c = 1:num_cell2D

        cell_index = py.int(c - 1);
        local_space_data = polydim.pde_tools.local_space_pcc_2_d.create_local_space(geometry_utilities_config.tolerance1_d, geometry_utilities_config.tolerance2_d, mesh_geometric_data, cell_index, reference_element_data);
        basis_functions_values = double(polydim.pde_tools.local_space_pcc_2_d.basis_functions_values(reference_element_data, local_space_data, polydim.vem.pcc.ProjectionTypes.pi0k));
        basis_functions_derivative_values = polydim.pde_tools.local_space_pcc_2_d.basis_functions_derivative_values(reference_element_data, local_space_data);
        basis_functions_derivative_values = {double(basis_functions_derivative_values{py.int(0)}), double(basis_functions_derivative_values{py.int(1)})};
        cell2_d_internal_quadrature = polydim.pde_tools.local_space_pcc_2_d.internal_quadrature(reference_element_data, local_space_data);

        weights = double(cell2_d_internal_quadrature.weights);
        points = double(cell2_d_internal_quadrature.points);
        exact_solution_values = exact_solution(points)';
        exact_derivative_solution_values = exact_derivative_solution(points)';

        global_do_fs = do_fs_data.cells_global_do_fs{cell_dimension}{cell_index};
        local_solution_do_fs = zeros(length(global_do_fs), 1);


        local_dof_i = int64(do_fs_data_indices.cells_do_fs_local_index{cell_index}) + 1;
        local_strong_i = int64(do_fs_data_indices.cells_strongs_local_index{cell_index}) + 1;

        global_dof_i = int64(do_fs_data_indices.cells_do_fs_global_index{cell_index}) + 1;
        global_strong_i = int64(do_fs_data_indices.cells_strongs_global_index{cell_index}) + 1;

        local_solution_do_fs(local_dof_i) = solution(global_dof_i);
        local_solution_do_fs(local_strong_i) = solution_dirichlet(global_strong_i);

        local_error_l2 = (basis_functions_values * local_solution_do_fs - exact_solution_values).^2;
        local_norm_l2 = (basis_functions_values * local_solution_do_fs).^2;

        cell2_ds_error_l2(c) = weights * local_error_l2;
        cell2_ds_norm_l2(c) = weights * local_norm_l2;

        local_error_h1 = ((basis_functions_derivative_values{1} * local_solution_do_fs - exact_derivative_solution_values{1}').^2 ...
            + (basis_functions_derivative_values{2} * local_solution_do_fs - exact_derivative_solution_values{2}').^2);

        local_norm_h1 = ((basis_functions_derivative_values{1} * local_solution_do_fs).^2 + ...
            (basis_functions_derivative_values{2} * local_solution_do_fs).^2);

        cell2_ds_error_h1(c) = sum(weights * local_error_h1);
        cell2_ds_norm_h1(c) = sum(weights * local_norm_h1);

        if mesh_geometric_data.cell2_ds_diameters{cell_index} > mesh_size
            mesh_size = mesh_geometric_data.cell2_ds_diameters{cell_index};
        end
    end


    error_l2 = sqrt(sum(cell2_ds_error_l2));
    norm_l2 = sqrt(sum(cell2_ds_norm_l2));
    error_h1 = sqrt(sum(cell2_ds_error_h1));
    norm_h1 = sqrt(sum(cell2_ds_norm_h1));

    errors = [errors; table(num_cell2D, num_dofs, num_strongs, mesh_size, error_l2, error_h1, norm_l2, norm_h1, residual_norm)];

    %% Export Solution
    cell0_ds_numeric = zeros(numCell0D, 1);
    cell0_ds_exact = exact_solution(cell0_ds_coordinates)';

    for p = 1:numCell0D

        local_do_fs = do_fs_data.cells_do_fs{py.int(0)}{py.int(p - 1)};

        local_dof_i = local_do_fs{py.int(0)};

        if local_dof_i.type == polydim.pde_tools.do_fs.DOFsManager.DOFsData.DOF.Types.strong
            cell0_ds_numeric(p) = solution_dirichlet(int64(local_dof_i.global_index) + 1);
        elseif local_dof_i.type == polydim.pde_tools.do_fs.DOFsManager.DOFsData.DOF.Types.dof
            cell0_ds_numeric(p) = solution(int64(local_dof_i.global_index) + 1);
        else
            error("Unknown DOF Type");
        end
    end

    vtk_utilities.export_solution_2(strcat([export_path, '/Solution_', int2str(int64(method_type.value)), '_', int2str(int64(method_order)), '_', int2str(ref)]), ...
        mesh, py.numpy.squeeze(py.numpy.array(cell0_ds_numeric, pyargs('dtype', py.numpy.float64))), py.numpy.squeeze(py.numpy.array(cell0_ds_exact, pyargs('dtype', py.numpy.float64))), ...
        py.numpy.squeeze(py.numpy.array(cell2_ds_error_l2, pyargs('dtype', py.numpy.float64))), ...
        py.numpy.squeeze(py.numpy.array(cell2_ds_error_h1, pyargs('dtype', py.numpy.float64))));

end

format short e

disp(errors)
writetable(errors,strcat([export_path, '/Errors_', int2str(int64(method_type.value)), '_', int2str(int64(method_order)), '_', int2str(ref), '.csv']),'Delimiter',';'); 

alpha = polyfit(log(errors.mesh_size), log(errors.error_l2), 1);
beta = polyfit(log(errors.mesh_size), log(errors.error_h1), 1);

disp('Convergence rates:')
disp([alpha(1) beta(1)])

%profile viewer

rmpath('Elliptic_PCC_2D/')