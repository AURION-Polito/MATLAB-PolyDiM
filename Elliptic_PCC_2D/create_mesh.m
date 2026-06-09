function [mesh, mesh_data, mesh_geometric_data] = create_mesh(gedim, polydim, geometry_utilities, mesh_utilities, mesh_type, mesh_max_relative_area, import_path, pde_domain)

mesh_data = gedim.MeshMatrices();
mesh = gedim.MeshMatricesDAO(mesh_data);

if polydim.pde_tools.mesh.pde_mesh_utilities.MeshGenerator_Types_2D.triangular == mesh_type || ...
        polydim.pde_tools.mesh.pde_mesh_utilities.MeshGenerator_Types_2D.minimal == mesh_type || ...
        polydim.pde_tools.mesh.pde_mesh_utilities.MeshGenerator_Types_2D.polygonal == mesh_type || ...
        polydim.pde_tools.mesh.pde_mesh_utilities.MeshGenerator_Types_2D.squared == mesh_type || ...
        polydim.pde_tools.mesh.pde_mesh_utilities.MeshGenerator_Types_2D.random_distorted == mesh_type
    polydim.pde_tools.mesh.pde_mesh_utilities.create_mesh_2_d(geometry_utilities, mesh_utilities, mesh_type, pde_domain, mesh_max_relative_area, mesh)
elseif polydim.pde_tools.mesh.pde_mesh_utilities.MeshGenerator_Types_2D.csv_importer == mesh_type || ...
        polydim.pde_tools.mesh.pde_mesh_utilities.MeshGenerator_Types_2D.off_importer == mesh_type
    polydim.pde_tools.mesh.pde_mesh_utilities.import_mesh_2_d(geometry_utilities, mesh_utilities, mesh_type, import_path, mesh)
else
    error("MeshGenerator " + str(mesh_type) + " not supported")
end


mesh_geometric_data_config = gedim.MeshUtilities.MeshGeometricData2DConfig(true, true, true, true, true, true, true, true, true, true, true, true, true);
mesh_geometric_data = polydim.pde_tools.mesh.pde_mesh_utilities.compute_mesh_2_d_geometry_data(geometry_utilities, mesh_utilities, mesh, mesh_geometric_data_config);

end