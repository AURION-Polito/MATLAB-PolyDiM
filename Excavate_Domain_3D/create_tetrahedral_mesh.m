function [mesh,mesh_data] = create_tetrahedral_mesh(gedim, domain_data, internal_points, domain_volume, mesh_max_volume, MAX_UINT)

mesh_data = gedim.MeshMatrices();
mesh = gedim.MeshMatricesDAO(mesh_data);

if isscalar(domain_data)
    gedim.external.mesh_utilities.create_tetrahedral_mesh(domain_data{1}.vertices,...
        domain_data{1}.edges,...
        domain_data{1}.faces,...
        max_volume, mesh);
else
    num_vertices = 0;
    num_faces = 0;
    for d = 1:length(domain_data)
        num_vertices = num_vertices + domain_data{d}.vertices.shape(2).int64;
        num_faces = num_faces + length(domain_data{d}.faces);
    end

    points = zeros(3, num_vertices);
    facets = py.list();
    offset_vertices = 1;
    offset_faces = 1;
    for d = 1:length(domain_data)

        num_local_vertices = domain_data{d}.vertices.shape(2).int64;
        points(:, offset_vertices:offset_vertices+num_local_vertices-1) = domain_data{d}.vertices;

        for f = 1:length(domain_data{d}.faces)
            local_facet = py.list();
            for v = 1:domain_data{d}.faces{f}.shape(2).int64
                face_vertices = int64(domain_data{d}.faces{f});
                local_facet.append(py.int(face_vertices(1, v) + offset_vertices - 1));
            end
            facets.append(local_facet);
            offset_faces = offset_faces + 1;
        end

        offset_vertices = offset_vertices + num_local_vertices;
    end

    holes_points = py.list();
    for h = 2:length(internal_points)
        holes_points.append(internal_points{h});
    end

    region = gedim.TetgenInterface.Region();
    region.centroid = internal_points{1};
    region.id = py.int(1);
    region.max_volume = mesh_max_volume * domain_volume;
    regions = {region};

    tetgen_interface = gedim.TetgenInterface();
    points = py.numpy.array(points, pyargs('dtype', py.numpy.float64));

    tetgen_interface.create_mesh(points, facets, regions, holes_points, mesh, 'QpqfezAa');

    for v = 1:int64(mesh.cell0_d_total_number())
        if mesh.cell0_d_marker(v-1) == MAX_UINT
            mesh.cell0_d_set_marker(v-1, num_faces + num_vertices + 1);
        end
    end

    for e = 1:int64(mesh.cell1_d_total_number())
        if mesh.cell1_d_marker(e-1) == MAX_UINT
            mesh.cell1_d_set_marker(e-1, num_faces + num_vertices + 1);

        end
    end
end
end