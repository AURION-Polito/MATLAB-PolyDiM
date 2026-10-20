%%
clc
clear
close all

currentPath = pwd;
file_path = fileparts(mfilename('fullpath'));

addpath(file_path + "/Excavate_Domain_3D/")

setup_pypolydim

%% Arguments

% Export Path
[~, ~, ~] = mkdir(currentPath + "/Export");
export_path = currentPath + "/Export/Excavate_Domain_3D";
[~, ~, ~] = mkdir(export_path);

export_vtk = true;
export_medit = true;
% Mesh Generator Type: 1 - PureTetrahedral; 2 - PureCubic;
mesh_type = 2;
% Mesh 3D maximum relative cell volume for tetrahedral mesh
max_volume = 0.1;
% Mesh 3D number of cells in each direction for parallelpipedal mesh. It
% must be a multiple of 4 for the tested case
num_cells_height = 12;
num_cells_base = 4;
num_cells_width = 16;

% Add internal points to edge to simulate edge dofs.
num_internal_points = 1;
% Internal Points Type:; 1 - Equispaced; 2 - Lobatto (Default: 1)
internal_points_type = 2;

% Geometric Tolerances
geometry_utilities_config = gedim.GeometryUtilitiesConfig();
geometry_utilities_config.tolerance1_d = 1.0e-12;
geometry_utilities_config.tolerance2_d = 1.0e-14;
geometry_utilities_config.tolerance3_d = 1.0e-15;
geometry_utilities = gedim.GeometryUtilities(geometry_utilities_config);


fprintf('\n');
fprintf('=========================================\n');
fprintf('      SIMULATION SETTINGS SUMMARY\n');
fprintf('=========================================\n');
fprintf(' Mesh Type                : %s\n', string(mesh_type));
fprintf(' Number of Internal Points: %d\n', num_internal_points);
fprintf(' Internal Point Type      : %s\n', string(internal_points_type));
fprintf('=========================================\n\n');

%% Domain Data
test_square_domain_with_one_hole

%% Generate Mesh
if mesh_type == 1
    [mesh, mesh_data] = create_tetrahedral_mesh(gedim, domain_data, internal_points, domain_volume, max_volume, MAX_UINT);
elseif mesh_type == 2
    [mesh, mesh_data] = create_parallelepid_mesh(gedim, geometry_utilities, domain_data, num_cells_base, num_cells_height, num_cells_width);
end

%% Add Hanging Nodes
add_edge_internal_points(gedim, num_internal_points, internal_points_type, mesh);

%% Export Mesh
if export_vtk
    export_mesh_path = export_path + "/MeshVtu";
    [~, ~, ~] = mkdir(export_mesh_path);
    ExportVTKUtilities.export_mesh(export_mesh_path, mesh)
    disp("Export VTU mesh in: " + export_mesh_path)
end

if export_medit
    export_mesh_path = export_path + "/MeshMEDIT";
    [~, ~, ~] = mkdir(export_mesh_path);
    gedim.MeshUtilities.compute_cell0_d_cell1_d_neighbours(mesh);
    gedim.MeshUtilities.export_mesh_to_medit(mesh, export_mesh_path, "Domain_Mesh", py.bool(0), py.int(num_internal_points + 1));
    disp("Export VTU mesh in: " + export_mesh_path)
end

rmpath(file_path + "/Excavate_Domain_3D/")