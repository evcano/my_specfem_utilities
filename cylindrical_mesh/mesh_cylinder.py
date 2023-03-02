#!/usr/bin/env python3

import fnmatch
import numpy as np
import pygmsh
import yaml
import sys


def mesh_cylinder(par):
    # READ PARAMETERS
    # ---------------
    mesh_name = par['MESH_NAME']
    x0 = par['X']
    y0 = par['Y']
    square_side = par['SQUARE_SIDE']
    square_curv = par['SQUARE_CURVATURE']
    circle_radius = par['CIRCLE_RADIUS']
    layers_depth = par['LAYERS_DEPTH']
    nex_theta = par['NEL_ANGULAR']
    nex_r = par['NEL_RADIAL']
    nex_layers = par['NEL_LAYER']

    # VERIFY PARAMETERS
    # -----------------
    if len(layers_depth)-1 != len(nex_layers):
        print('Please define the number of elements in each layer\n')
        sys.exit()

    if len(circle_radius) != len(nex_r):
        print('Please define the number of elements in each circle\n')
        sys.exit()

    if (nex_theta % 4) != 0:
        print('NEL_ANGULAR must be a multiple of 4\n')
        sys.exit()

    # SET SOME VARIABLES
    # ------------------
    mesh_algorithm = 8
    mesh_order = 1
    mesh_type = 'Progression'
    coeff = 1.0
    arrangement = 'Left'

    ncircles = len(circle_radius)
    nlayers = len(layers_depth)

    nex_theta = int(nex_theta / 4)
    nex_theta += 1

    for i in range(0, ncircles):
        nex_r[i] += 1

    for i in range(0, nlayers-1):
        nex_layers[i] += 1

    point = {}
    line = {}
    surface = {}
    volume = {}

    # INITIALIZE MESH
    # ---------------
    G = pygmsh.geo.Geometry()
    MO = G.__enter__()

    # DEFINE LAYERS
    # -------------
    for z, depth in enumerate(layers_depth):
        # POINTS
        # ------
        # center
        center_point_coor = np.array([x0, y0, depth])
        ep = MO.add_point(center_point_coor)
        point[f'center_{z}'] = ep

        # square
        a = square_side / 2.0
        b = a * square_curv

        point_coor = [center_point_coor + np.array([-a, a, 0.0]),
                      center_point_coor + np.array([0.0, b, 0.0]),
                      center_point_coor + np.array([a, a, 0.0]),
                      center_point_coor + np.array([b, 0.0, 0.0]),
                      center_point_coor + np.array([a, -a, 0.0]),
                      center_point_coor + np.array([0.0, -b, 0.0]),
                      center_point_coor + np.array([-a, -a, 0.0]),
                      center_point_coor + np.array([-b, 0.0, 0.0])]

        point[f'square_{z}'] = []

        for coor in point_coor:
            ep = MO.add_point(coor)
            point[f'square_{z}'].append(ep)

        # circles
        for r, radius in enumerate(circle_radius):
            a = radius / np.sqrt(2)

            point_coor = [center_point_coor + np.array([-a, a, 0.0]),
                          center_point_coor + np.array([a, a, 0.0]),
                          center_point_coor + np.array([a, -a, 0.0]),
                          center_point_coor + np.array([-a, -a, 0.0])]

            point[f'circle_{z}_{r}'] = []

            for coor in point_coor:
                ep = MO.add_point(coor)
                point[f'circle_{z}_{r}'].append(ep)

        # LINES
        # -----
        # square sides
        line[f'square_{z}'] = []

        for i in range(0, 7, 2):
            j = i + 2 if i != 6 else 0

            el = MO.add_spline([point[f'square_{z}'][i],
                                point[f'square_{z}'][i+1],
                                point[f'square_{z}'][j]])

            line[f'square_{z}'].append(el)

        # circles arcs
        for r in range(0, ncircles):
            line[f'circle_{z}_{r}'] = []

            for i in range(0, 4):
                j = i + 1 if i != 3 else 0

                el = MO.add_circle_arc(point[f'circle_{z}_{r}'][i],
                                       point[f'center_{z}'],
                                       point[f'circle_{z}_{r}'][j])

                line[f'circle_{z}_{r}'].append(el)

        # join first circle with square
        line[f'joint_{z}_0'] = []

        for j, i in enumerate(range(0, 7, 2)):
            el = MO.add_line(point[f'square_{z}'][i], point[f'circle_{z}_0'][j])
            line[f'joint_{z}_0'].append(el)

        # join other circles
        for r in range(1, ncircles):
            line[f'joint_{z}_{r}'] = []

            for i in range(0, 4):
                el = MO.add_line(point[f'circle_{z}_{r-1}'][i],
                                 point[f'circle_{z}_{r}'][i])

                line[f'joint_{z}_{r}'].append(el)

        # SURFACES
        # --------
        # square
        ello = MO.add_curve_loop(line[f'square_{z}'])
        es = MO.add_surface(ello)
        surface[f'square_{z}'] = es

        # circles
        for r in range(0, ncircles):
            surface[f'circle_{z}_{r}'] = []

            for i in range(0, 4):
                j = i + 1 if i != 3 else 0

                if r == 0:
                    ello = MO.add_curve_loop([line[f'circle_{z}_{r}'][i],
                                             -line[f'joint_{z}_{r}'][j],
                                             -line[f'square_{z}'][i],
                                             line[f'joint_{z}_{r}'][i]])
                else:
                    ello = MO.add_curve_loop([line[f'circle_{z}_{r}'][i],
                                             -line[f'joint_{z}_{r}'][j],
                                             -line[f'circle_{z}_{r-1}'][i],
                                             line[f'joint_{z}_{r}'][i]])

                es = MO.add_surface(ello)
                surface[f'circle_{z}_{r}'].append(es)

    # DEFINE SIDES
    # ------------
    for z in range(0, nlayers-1):
        # LINES
        # -----
        # square
        line[f'side_square_{z}'] = []

        for i in range(0, 7, 2):
            el = MO.add_line(point[f'square_{z}'][i],
                             point[f'square_{z+1}'][i])

            line[f'side_square_{z}'].append(el)

        # circle
        for r in range(0, ncircles):
            line[f'side_circle_{z}_{r}'] = []

            for i in range(0, 4):
                el = MO.add_line(point[f'circle_{z}_{r}'][i],
                                 point[f'circle_{z+1}_{r}'][i])

                line[f'side_circle_{z}_{r}'].append(el)

        # SURFACES
        # --------
        # square
        surface[f'side_square_{z}'] = []

        for i in range(0, 4):
            j = i + 1 if i != 3 else 0

            ello = MO.add_curve_loop([line[f'square_{z}'][i],
                                     line[f'side_square_{z}'][j],
                                     -line[f'square_{z+1}'][i],
                                     -line[f'side_square_{z}'][i]])

            es = MO.add_surface(ello)
            surface[f'side_square_{z}'].append(es)

        # circles
        for r in range(0, ncircles):
            surface[f'side_circle_{z}_{r}'] = []

            for i in range(0, 4):
                j = i + 1 if i != 3 else 0

                ello = MO.add_curve_loop([line[f'circle_{z}_{r}'][i],
                                          line[f'side_circle_{z}_{r}'][j],
                                          -line[f'circle_{z+1}_{r}'][i],
                                          -line[f'side_circle_{z}_{r}'][i]])

                es = MO.add_surface(ello)
                surface[f'side_circle_{z}_{r}'].append(es)

        # joints
        for r in range(0, ncircles):
            surface[f'side_joint_{z}_{r}'] = []

            for i in range(0, 4):
                if r == 0:
                    ello = MO.add_curve_loop([-line[f'joint_{z}_{r}'][i],
                                              line[f'side_square_{z}'][i],
                                              line[f'joint_{z+1}_{r}'][i],
                                              -line[f'side_circle_{z}_{r}'][i]])
                else:
                    ello = MO.add_curve_loop([-line[f'joint_{z}_{r}'][i],
                                              line[f'side_circle_{z}_{r-1}'][i],
                                              line[f'joint_{z+1}_{r}'][i],
                                              -line[f'side_circle_{z}_{r}'][i]])

                es = MO.add_surface(ello)
                surface[f'side_joint_{z}_{r}'].append(es)

        # VOLUMES
        # --------
        # square
        eslo = MO.add_surface_loop([surface[f'square_{z}'],
                                    surface[f'side_square_{z}'][0],
                                    surface[f'side_square_{z}'][1],
                                    surface[f'side_square_{z}'][2],
                                    surface[f'side_square_{z}'][3],
                                    surface[f'square_{z+1}']])
        ev = MO.add_volume(eslo)
        volume[f'square_{z}'] = ev

        # circles
        for r in range(0, ncircles):
            volume[f'circle_{z}_{r}'] = []

            for i in range(0, 4):
                j = i + 1 if i != 3 else 0

                if r == 0:
                    eslo = MO.add_surface_loop([surface[f'circle_{z}_{r}'][i],
                                                surface[f'side_square_{z}'][i],
                                                surface[f'circle_{z+1}_{r}'][i],
                                                surface[f'side_circle_{z}_{r}'][i],
                                                surface[f'side_joint_{z}_{r}'][i],
                                                surface[f'side_joint_{z}_{r}'][j]])
                else:
                    eslo = MO.add_surface_loop([surface[f'circle_{z}_{r}'][i],
                                                surface[f'side_circle_{z}_{r-1}'][i],
                                                surface[f'circle_{z+1}_{r}'][i],
                                                surface[f'side_circle_{z}_{r}'][i],
                                                surface[f'side_joint_{z}_{r}'][i],
                                                surface[f'side_joint_{z}_{r}'][j]])

                ev = MO.add_volume(eslo)
                volume[f'circle_{z}_{r}'].append(ev)

    # TRANSFINITE MESH
    # ----------------
    # lines
    for key, lines in line.items():
        if fnmatch.fnmatch(key, 'joint*'):
            r = int(key[-1])
            num_nodes = nex_r[r]
        elif fnmatch.fnmatch(key, 'side_square*'):
            z = int(key[-1])
            num_nodes = nex_layers[z]
        elif fnmatch.fnmatch(key, 'side_*_*_*'):
            z = int(key[-3])
            num_nodes = nex_layers[z]
        else:
            num_nodes = nex_theta

        for el in lines:
            MO.set_transfinite_curve(curve=el, num_nodes=num_nodes,
                                     mesh_type=mesh_type, coeff=coeff)

    # surfaces
    for surfaces in surface.values():
        if not isinstance(surfaces, list):
            surfaces = [surfaces]

        for es in surfaces:
            MO.set_transfinite_surface(surface=es, arrangement=arrangement,
                                       corner_pts=[])

        MO.set_recombined_surfaces(surfaces)

    # volumes
    for volumes in volume.values():
        if not isinstance(volumes, list):
            volumes = [volumes]

        for ev in volumes:
            MO.set_transfinite_volume(volume=ev, corner_pts=[])

    # DEFINE PHYSICAL VOLUMES
    # -------------------------
    for z in range(0, nlayers -1):
        label = f'm{z+1}'
        layer = []
        layer.append(volume[f'square_{z}'])

        for r in range(0, ncircles):
            layer.extend(volume[f'circle_{z}_{r}'])

        MO.add_physical(layer, label=label)

    # DEFINE PHYSICAL SURFACES
    # -------------------------
    top = []
    bottom = []

    top.append(surface['square_0'])
    bottom.append(surface[f'square_{nlayers-1}'])

    for r in range(0, ncircles):
        top.extend(surface[f'circle_0_{r}'])
        bottom.extend(surface[f'circle_{nlayers-1}_{r}'])

    MO.add_physical(top, label='top')
    MO.add_physical(bottom, label='bottom')

    xmin = []
    xmax = []
    ymin = []
    ymax = []

    for z in range(0, nlayers - 1):
        xmin.append(surface[f'side_circle_{z}_{ncircles-1}'][3])
        xmax.append(surface[f'side_circle_{z}_{ncircles-1}'][1])

        ymin.append(surface[f'side_circle_{z}_{ncircles-1}'][2])
        ymax.append(surface[f'side_circle_{z}_{ncircles-1}'][0])

    MO.add_physical(xmin, label='xmin')
    MO.add_physical(xmax, label='xmax')
    MO.add_physical(ymin, label='ymin')
    MO.add_physical(ymax, label='ymax')

    # GENERATE AND SAVE MESH
    # ----------------------
    ME = MO.generate_mesh(algorithm=mesh_algorithm,
                          order=mesh_order,
                          verbose=True)

    pygmsh.write(mesh_name)

    return


if __name__ == '__main__':
    if len(sys.argv) != 2:
        print('usage: ./cylinder_gmsh.py file\n')
        sys.exit(1)
    else:
        parfile_name = sys.argv[1]

    try:
        with open(parfile_name, 'r') as _file:
            par = yaml.load(_file)
    except IOError:
        print('IOError: {} not found.'.format(parfile_name))

    mesh_cylinder(par)
