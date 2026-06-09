function global_rhs = compute_weak_term(polydim, gedim, cell3_d_index, mesh, mesh_geometric_data, mesh_do_fs_info, do_fs_data, reference_element_data, local_space_data, weak_boundary_condition, global_rhs)

%% ------------------------------------------------------------------------
% Assembles the weak (Neumann) boundary contributions associated with the
% current 2D cell and adds them to the global right-hand side vector.
%
% The boundary integral is evaluated by Gaussian quadrature on each edge
% subject to a weak boundary condition. The resulting contributions are
% projected onto the local basis functions and accumulated in the global
% right-hand side.
%
% INPUTS:
%   polydim                : PolyDim library.
%
%   gedim                  : GeDiM library.
%
%   cell3_d_index          : Global index of the current 3D cell.
%
%   mesh                   : Mesh object containing the geometric and
%                            topological information of the computational
%                            domain.
%
%   mesh_geometric_data    : Structure containing the geometric quantities
%                            associated with the mesh, such as vertices,
%                            edge tangents, lengths, and orientations.
%
%   mesh_do_fs_info        : Structure containing the boundary
%                            classification and DOF-related information
%                            for the mesh entities.
%
%   do_fs_data             : Structure storing the degrees of freedom
%                            associated with mesh entities.
%
%   reference_element_data : Data defining the reference finite element,
%                            including the polynomial order and local
%                            geometry.
%
%   local_space_data       : Data describing the local approximation space
%                            and its associated basis functions.
%
%   weak_boundary_condition: Function handle that evaluates the prescribed
%                            Neumann boundary condition. It takes as input
%                            a boundary marker and the coordinates of the
%                            evaluation points and returns the
%                            corresponding boundary values.
%
%   global_rhs             : Global right-hand side vector to which the
%                            weak boundary contributions are added.
%
% OUTPUTS:
%   The function updates the vector `global_rhs` by adding the
%   contributions arising from the weak (Neumann) boundary conditions.
%% ------------------------------------------------------------------------



num_faces = int64(mesh.cell3_d_number_faces(cell3_d_index));
quadrature_point_offset = py.int(0);
for f = 1:num_faces

    cell2_d_index = mesh.cell3_d_face(cell3_d_index, py.int(f - 1));
    face_quadrature = polydim.pde_tools.local_space_pcc_3_d.face_quadrature(reference_element_data, local_space_data, py.int(f - 1), quadrature_point_offset);
    boundary_info = mesh_do_fs_info.cells_boundary_info{py.int(2)}{cell2_d_index};

    if boundary_info.type ~= polydim.pde_tools.do_fs.DOFsManager.MeshDOFsInfo.BoundaryInfo.BoundaryTypes.weak
        continue
    end
    
    weak_quadrature_points = double(face_quadrature.points);
    weak_quadrature_weights = double(face_quadrature.weights);
    neumann_values = weak_boundary_condition(int64(boundary_info.marker), weak_quadrature_points)';
    weak_basis_function_values = double(polydim.pde_tools.local_space_pcc_3_d.basis_functions_values_on_face(py.int(f - 1), reference_element_data, local_space_data, face_quadrature.points));

    % compute values of Neumann
    neumann_contributions = weak_basis_function_values' * diag(weak_quadrature_weights) * neumann_values;
    
    offset_pi0_km1 = 0;
    for p = 1:int64(mesh.cell2_d_number_vertices(cell2_d_index))
        cell0_d_index = mesh.cell2_d_vertex(cell2_d_index, py.int(p - 1));
        local_do_fs = do_fs_data.cells_do_fs{py.int(0)}{cell0_d_index};

        local_dof_i = local_do_fs{py.int(0)};

        if local_dof_i.type == polydim.pde_tools.do_fs.DOFsManager.DOFsData.DOF.Types.strong
        elseif local_dof_i.type == polydim.pde_tools.do_fs.DOFsManager.DOFsData.DOF.Types.dof
            global_rhs(int64(local_dof_i.global_index) + 1) = global_rhs(int64(local_dof_i.global_index) + 1) + neumann_contributions(offset_pi0_km1 + 1);
        else
            error("Unknown DOF Type")
        end
        
        offset_pi0_km1 = offset_pi0_km1 + 1;
    end

    for p = 1:int64(mesh.cell2_d_number_edges(cell2_d_index))
        cell1_d_index = mesh.cell2_d_edge(cell2_d_index, py.int(p - 1));
        local_do_fs = do_fs_data.cells_do_fs{py.int(1)}{cell1_d_index};
        
        for loc_i = 1:length(local_do_fs)
            local_dof_i = local_do_fs{py.int(loc_i - 1)};
    
            if local_dof_i.type == polydim.pde_tools.do_fs.DOFsManager.DOFsData.DOF.Types.strong
            elseif local_dof_i.type == polydim.pde_tools.do_fs.DOFsManager.DOFsData.DOF.Types.dof
                global_rhs(int64(local_dof_i.global_index) + 1) = global_rhs(int64(local_dof_i.global_index) + 1) + neumann_contributions(loc_i + offset_pi0_km1);
            else
                error("Unknown DOF Type")
            end
        end

        offset_pi0_km1 = offset_pi0_km1 + length(local_do_fs);
    end

    local_do_fs = do_fs_data.cells_do_fs{py.int(2)}{cell2_d_index};

    for loc_i = 1:length(local_do_fs)
        local_dof_i = local_do_fs{py.int(loc_i - 1)};

        if local_dof_i.type == polydim.pde_tools.do_fs.DOFsManager.DOFsData.DOF.Types.strong
        elseif local_dof_i.type == polydim.pde_tools.do_fs.DOFsManager.DOFsData.DOF.Types.dof
            global_rhs(int64(local_dof_i.global_index) + 1) = global_rhs(int64(local_dof_i.global_index) + 1) + neumann_contributions(loc_i + offset_pi0_km1);
        else
            error("Unknown DOF Type")
        end
    end
end

end