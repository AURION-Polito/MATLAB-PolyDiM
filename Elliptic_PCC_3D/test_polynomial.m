pde_domain = polydim.pde_tools.mesh.pde_mesh_utilities.PDE_Domain_3D();

% set verticse
pde_domain.vertices = [0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0;
                        0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0;
                        0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0];

% set edges
pde_domain.edges = int64([0, 1, 2, 3, 4, 5, 6, 7, 0, 1, 2, 3;
                          1, 2, 3, 0, 5, 6, 7, 4, 4, 5, 6, 7]);

% set faces
pde_domain.faces = { int64([0, 1, 2, 3; 0, 1, 2, 3]), int64([4, 5, 6, 7; 4, 5, 6, 7]), ...
                     int64([0, 3, 7, 4; 3, 11, 7, 8]), int64([1, 2, 6, 5; 1, 10, 5, 9]), ...
                     int64([0, 1, 5, 4; 0, 9, 4, 8]), int64([3, 2, 6, 7;2, 10, 6, 11])};

pde_domain.volume = 1.0;
pde_domain.shape_type = polydim.pde_tools.mesh.pde_mesh_utilities.PDE_Domain_3D.Domain_Shape_Types.parallelepiped;

info_internal = polydim.pde_tools.do_fs.DOFsManager.MeshDOFsInfo.BoundaryInfo(polydim.pde_tools.do_fs.DOFsManager.MeshDOFsInfo.BoundaryInfo.BoundaryTypes.none);
info_internal.marker = py.int(0);

info_neumann_bottom = polydim.pde_tools.do_fs.DOFsManager.MeshDOFsInfo.BoundaryInfo(polydim.pde_tools.do_fs.DOFsManager.MeshDOFsInfo.BoundaryInfo.BoundaryTypes.weak);
info_neumann_bottom.marker = py.int(2);

info_dirichlet = polydim.pde_tools.do_fs.DOFsManager.MeshDOFsInfo.BoundaryInfo(polydim.pde_tools.do_fs.DOFsManager.MeshDOFsInfo.BoundaryInfo.BoundaryTypes.strong);
info_dirichlet.marker = py.int(1);

boundary_info = py.dict();

boundary_info{py.int(0)} = info_internal;

boundary_info{py.int(1)} = info_dirichlet;
boundary_info{py.int(2)} = info_dirichlet;
boundary_info{py.int(3)} = info_dirichlet;
boundary_info{py.int(4)} = info_dirichlet;
boundary_info{py.int(5)} = info_dirichlet;
boundary_info{py.int(6)} = info_dirichlet;
boundary_info{py.int(7)} = info_dirichlet;
boundary_info{py.int(8)} = info_dirichlet;

boundary_info{py.int(9)} = info_dirichlet;
boundary_info{py.int(10)} = info_dirichlet;
boundary_info{py.int(11)} = info_dirichlet;
boundary_info{py.int(12)} = info_dirichlet;
boundary_info{py.int(13)} = info_dirichlet;
boundary_info{py.int(14)} = info_dirichlet;
boundary_info{py.int(15)} = info_dirichlet;
boundary_info{py.int(16)} = info_dirichlet;
boundary_info{py.int(17)} = info_dirichlet;
boundary_info{py.int(18)} = info_dirichlet;
boundary_info{py.int(19)} = info_dirichlet;
boundary_info{py.int(20)} = info_dirichlet;

boundary_info{py.int(21)} = info_neumann_bottom;
boundary_info{py.int(22)} = info_dirichlet;
boundary_info{py.int(23)} = info_dirichlet;
boundary_info{py.int(24)} = info_dirichlet;
boundary_info{py.int(25)} = info_dirichlet;
boundary_info{py.int(26)} = info_dirichlet;

diffusion_term = @(points) ones(1, size(points, 2));
source_term = @(points) 128.0 * (points(1, :) .* (1.0 - points(1, :)) .* points(3, :) .* (1.0 - points(3, :)) ...
    + points(2, :) .* (1.0 - points(2, :)) .* points(1, :) .* (1.0 - points(1, :)) ...
    + points(2, :) .* (1.0 - points(2, :)) .* points(3, :) .* (1.0 - points(3, :)));

strong_boundary_condition = @(marker, points) (marker == 1) * (64.0 * points(1, :) .* (1.0 - points(1, :)) .* points(2, :) .* (1.0 - points(2, :)) .* points(3, :) .* (1.0 - points(3, :)) + 1.1);

weak_boundary_condition = @(marker, points) -(marker == 2) * (64.0 * points(1, :) .* (1.0 - points(1, :)) .* (1.0 - 2.0 * points(3, :)) .* points(2, :) .* (1.0 - points(2, :)));

exact_solution = @(points) (64.0 * points(1, :) .* (1.0 - points(1, :)) .* points(2, :) .* (1.0 - points(2, :)) .* points(3, :) .* (1.0 - points(3, :)) + 1.1);

exact_derivative_solution = @(points) {64.0 * (1.0 - 2.0 * points(1, :)) .* points(2, :) .* (1.0 - points(2, :)) .* points(3, :) .* (1.0 - points(3, :));
    64.0 * points(1, :) .* (1.0 - points(1, :)) .* (1.0 - 2.0 * points(2, :)) .* points(3, :) .* (1.0 - points(3, :));
    64.0 * points(1, :) .* (1.0 - points(1, :)) .* (1.0 - 2.0 * points(3, :)) .* points(2, :) .* (1.0 - points(2, :))};