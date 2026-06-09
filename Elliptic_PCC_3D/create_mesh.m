function [mesh, mesh_data, mesh_geometric_data] = create_mesh(gedim, polydim, geometry_utilities, mesh_utilities, mesh_type, mesh_max_relative_volume, import_path, pde_domain)

mesh_data = gedim.MeshMatrices();
mesh = gedim.MeshMatricesDAO(mesh_data);

if polydim.pde_tools.mesh.pde_mesh_utilities.MeshGenerator_Types_3D.tetrahedral == mesh_type || ...
        polydim.pde_tools.mesh.pde_mesh_utilities.MeshGenerator_Types_3D.minimal == mesh_type || ...
        polydim.pde_tools.mesh.pde_mesh_utilities.MeshGenerator_Types_3D.polyhedral == mesh_type || ...
        polydim.pde_tools.mesh.pde_mesh_utilities.MeshGenerator_Types_3D.cubic == mesh_type 
    polydim.pde_tools.mesh.pde_mesh_utilities.create_mesh_3_d(geometry_utilities, mesh_utilities, mesh_type, pde_domain, mesh_max_relative_volume, mesh)
elseif polydim.pde_tools.mesh.pde_mesh_utilities.MeshGenerator_Types_3D.csv_importer == mesh_type || ...
        polydim.pde_tools.mesh.pde_mesh_utilities.MeshGenerator_Types_3D.vtk_importer == mesh_type || ...
        polydim.pde_tools.mesh.pde_mesh_utilities.MeshGenerator_Types_3D.off_importer == mesh_type
    polydim.pde_tools.mesh.pde_mesh_utilities.import_mesh_3_d(mesh_utilities, mesh_type, import_path, mesh)
else
    error("MeshGenerator " + str(mesh_type) + " not supported")
end

mesh_geometric_data = polydim.pde_tools.mesh.pde_mesh_utilities.compute_mesh_3_d_geometry_data(geometry_utilities, mesh);

end