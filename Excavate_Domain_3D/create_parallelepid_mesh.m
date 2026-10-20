function [mesh, mesh_data] = create_parallelepid_mesh(gedim, geometry_utilities, domain_data, num_cells_base, num_cells_height, num_cells_width)

if mod(num_cells_base, 4) ~= 0 || mod(num_cells_height, 4) ~= 0 || mod(num_cells_width, 4) ~= 0
    error('not valid number of cells in each dimension')
end

mesh_data = gedim.MeshMatrices();
mesh = gedim.MeshMatricesDAO(mesh_data);

vertices = double(domain_data{1}.vertices);
domain_origin = py.numpy.array(vertices(:, 1), pyargs('dtype', py.numpy.float64));
domain_base_tangent = py.numpy.array(vertices(:, 2) - domain_origin, pyargs('dtype', py.numpy.float64));
domain_width_tangent = py.numpy.array(vertices(:, 5) - domain_origin, pyargs('dtype', py.numpy.float64));
domain_heigth_tangent = py.numpy.array(vertices(:, 4) - domain_origin, pyargs('dtype', py.numpy.float64));

curvilinear_base = py.numpy.array(linspace(0, 1, num_cells_base + 1), pyargs('dtype', py.numpy.float64));
curvilinear_height = py.numpy.array(linspace(0, 1, num_cells_height + 1), pyargs('dtype', py.numpy.float64));
curvilinear_width = py.numpy.array(linspace(0, 1, num_cells_width + 1), pyargs('dtype', py.numpy.float64));

if isscalar(domain_data)
    gedim.MeshUtilities.create_parallelepiped_mesh(domain_origin,...
        domain_base_tangent, domain_heigth_tangent, domain_width_tangent,...
        curvilinear_base, curvilinear_height, curvilinear_width, mesh);
else
    full_mesh_data = gedim.MeshMatrices();
    full_mesh = gedim.MeshMatricesDAO(full_mesh_data);
    gedim.MeshUtilities.create_parallelepiped_mesh(domain_origin,...
        domain_base_tangent, domain_heigth_tangent, domain_width_tangent,...
        curvilinear_base, curvilinear_height, curvilinear_width, full_mesh);

     gedim.external.mesh_utilities.create_mesh_with_holes(geometry_utilities, domain_data, full_mesh_data, full_mesh, mesh);
end

end