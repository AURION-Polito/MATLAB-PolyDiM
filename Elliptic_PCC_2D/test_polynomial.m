pde_domain = polydim.pde_tools.mesh.pde_mesh_utilities.PDE_Domain_2D();

pde_domain.vertices = [0.0, 1.0, 1.0, 0.0;
    0.0, 0.0, 1.0, 1.0;
    0.0, 0.0, 0.0, 0.0];

pde_domain.area = 1.0;
pde_domain.shape_type = polydim.pde_tools.mesh.pde_mesh_utilities.PDE_Domain_2D.Domain_Shape_Types.parallelogram;

info_internal = polydim.pde_tools.do_fs.DOFsManager.MeshDOFsInfo.BoundaryInfo(polydim.pde_tools.do_fs.DOFsManager.MeshDOFsInfo.BoundaryInfo.BoundaryTypes.none);
info_internal.marker = py.int(0);

info_neumann_right = polydim.pde_tools.do_fs.DOFsManager.MeshDOFsInfo.BoundaryInfo(polydim.pde_tools.do_fs.DOFsManager.MeshDOFsInfo.BoundaryInfo.BoundaryTypes.weak);
info_neumann_right.marker = py.int(2);

info_dirichlet = polydim.pde_tools.do_fs.DOFsManager.MeshDOFsInfo.BoundaryInfo(polydim.pde_tools.do_fs.DOFsManager.MeshDOFsInfo.BoundaryInfo.BoundaryTypes.strong);
info_dirichlet.marker = py.int(1);

boundary_info = py.dict();

boundary_info{py.int(0)} = info_internal;
boundary_info{py.int(1)} = info_dirichlet;
boundary_info{py.int(2)} = info_dirichlet;
boundary_info{py.int(3)} = info_dirichlet;
boundary_info{py.int(4)} = info_dirichlet;

boundary_info{py.int(5)} = info_dirichlet;
boundary_info{py.int(6)} = info_neumann_right;
boundary_info{py.int(7)} = info_dirichlet;
boundary_info{py.int(8)} = info_dirichlet;

diffusion_term = @(points) ones(1, size(points, 2));
source_term = @(points) 32.0 * (points(1, :) .* (1.0 - points(1, :)) + points(2, :) .* (1.0 - points(2, :)));
strong_boundary_condition = @(marker, points) (marker == 1) * (16.0 * points(1, :) .* (1.0 - points(1, :)) .* points(2, :) .* (1.0 - points(2, :)) + 1.1);
weak_boundary_condition = @(marker, points) (marker == 2) * (16.0 * (1.0 - 2.0 * points(1, :)) .* points(2, :) .* (1.0 - points(2, :)));
exact_solution = @(points) (16.0 * points(1, :) .* (1.0 - points(1, :)) .* points(2, :) .* (1.0 - points(2, :)) + 1.1);
exact_derivative_solution = @(points) {16.0 * (1.0 - 2.0 * points(1, :)) .* points(2, :) .* (1.0 - points(2, :));
    16.0 * points(1, :) .* (1.0 - points(1, :)) .* (1.0 - 2.0 * points(2, :))};