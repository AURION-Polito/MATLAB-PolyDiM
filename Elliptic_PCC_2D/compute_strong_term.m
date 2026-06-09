function solution_dirichlet = compute_strong_term(polydim, cell2_d_index, mesh, mesh_do_fs_info, do_fs_data, reference_element_data, local_space_data, strong_boundary_condition, solution_dirichlet)

%% ------------------------------------------------------------------------
% Assemble strong (Dirichlet) boundary conditions for the current 2D cell.
%
% The prescribed values are stored in the vector solution_dirichlet by
% assigning the corresponding global DOFs associated with strong boundary
% conditions.
%
% INPUTS:
%   polydim                   : PolyDim library
%
%   cell2_d_index             : Global index of the current 2D cell.
%
%   mesh                      : Mesh object containing the geometric and
%                               topological information of the computational
%                               domain.
%
%   mesh_do_fs_info           : Structure containing the boundary
%                               classification and DOF-related information
%                               for the mesh entities.
%
%   do_fs_data                : Structure storing the degrees of freedom
%                               associated with mesh entities.
%
%   reference_element_data    : Data defining the reference element,
%                               including the local geometry used for DOF
%                               placement.
%
%   local_space_data          : Data describing the local approximation
%                               space and its associated basis functions.
%
%   strong_boundary_condition : Function handle that evaluates the
%                               prescribed Dirichlet boundary condition.
%                               It takes as input a boundary marker and the
%                               coordinates of the evaluation points and
%                               returns the corresponding boundary values.
%
%   solution_dirichlet        : Global vector where the values associated
%                               with strong boundary DOFs are assigned.
%
% OUTPUTS:
%   The function updates the vector `solution_dirichlet` by assigning the
%   values corresponding to strong boundary conditions.
%% ------------------------------------------------------------------------


% Assemble strong boundary condition on Cell0Ds
coordinates = double(mesh.cell2_d_vertices_coordinates(cell2_d_index));
for v = 0:mesh.cell2_d_number_vertices(cell2_d_index)-1

    cell0_d_index = mesh.cell2_d_vertex(cell2_d_index, py.int(v));
    boundary_info = mesh_do_fs_info.cells_boundary_info{py.int(0)}{cell0_d_index};

    if boundary_info.type ~= polydim.pde_tools.do_fs.DOFsManager.MeshDOFsInfo.BoundaryInfo.BoundaryTypes.strong
        continue
    end

    strong_boundary_values = strong_boundary_condition(int64(boundary_info.marker), coordinates(:, v + 1));

    local_dofs = do_fs_data.cells_do_fs{py.int(0)}{cell0_d_index};

    for loc_i=1:length(local_dofs)

        local_dof_i = local_dofs{py.int(loc_i - 1)};

        if local_dof_i.type == polydim.pde_tools.do_fs.DOFsManager.DOFsData.DOF.Types.strong
            solution_dirichlet(int64(local_dof_i.global_index) + 1) = strong_boundary_values(loc_i);
        elseif local_dof_i.type == polydim.pde_tools.do_fs.DOFsManager.DOFsData.DOF.Types.dof
            continue
        else
            error("Unknown DOF Type")
        end

    end
end

% Assemble strong boundary condition on Cell1Ds
for ed=0:mesh.cell2_d_number_edges(cell2_d_index)-1

    cell1_d_index = mesh.cell2_d_edge(cell2_d_index, py.int(ed));

    boundary_info = mesh_do_fs_info.cells_boundary_info{py.int(1)}{cell1_d_index};
    local_dofs = do_fs_data.cells_do_fs{py.int(1)}{cell1_d_index};

    if boundary_info.type ~= polydim.pde_tools.do_fs.DOFsManager.MeshDOFsInfo.BoundaryInfo.BoundaryTypes.strong || isempty(local_dofs) 
        continue
    end

    edge_do_fs_coordinates = double(polydim.pde_tools.local_space_pcc_2_d.edge_dofs_coordinates(reference_element_data, local_space_data, py.int(ed)));

    strong_boundary_values = strong_boundary_condition(int64(boundary_info.marker), edge_do_fs_coordinates);

    for loc_i = 1:length(local_dofs)

        local_dof_i = local_dofs{py.int(loc_i - 1)};

        if local_dof_i.type == polydim.pde_tools.do_fs.DOFsManager.DOFsData.DOF.Types.strong
            solution_dirichlet(local_dof_i.global_index + 1) = strong_boundary_values(loc_i);
        elseif local_dof_i.type == polydim.pde_tools.do_fs.DOFsManager.DOFsData.DOF.Types.dof
            continue
        else
            error("Unknown DOF Type")
        end
    end
end

end