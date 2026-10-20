function [] = add_edge_internal_points(gedim, num_internal_points, internal_points_type, mesh)

if num_internal_points == 0
    return;
end

hanging_per_edges = py.dict();

if internal_points_type == 1
    v = linspace(0, 1, num_internal_points+2);
    reference_points = py.numpy.array(v(2:end-1), pyargs('dtype', py.numpy.float64));
    for e = 0:int64(mesh.cell1_d_total_number())-1
        hanging_per_edges{e} = reference_points;
    end
elseif internal_points_type == 2
    reference_data = gedim.quadrature.Quadrature_GaussLobatto1D.fill_points_and_weights(py.int(2 * num_internal_points + 1));
    ref_pts = double(reference_data.points);
    
    reference_points = py.list();
    for i = 1:num_internal_points
        reference_points.append(ref_pts(1, i + 1));
    end

    for e = 0:int64(mesh.cell1_d_total_number())-1
       hanging_per_edges{e} = reference_points;
    end
else
    error("not valid type: (1) equispaced - (2) lobatto");
end

gedim.MeshUtilities.add_hanging_nodes_to_mesh(hanging_per_edges, mesh);

end