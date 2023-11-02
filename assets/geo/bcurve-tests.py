import numpy as np
import scipy.integrate as sci

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os, math, json
from itertools import chain

from geomdl import NURBS
from geomdl import helpers
from geomdl import utilities
from geomdl import operations
# Import Matplotlib visualization module

file_name = 'bcurve-tests.json'
# max order
pmax = 4
# knots
knots_maxp = [0,0,0,0,0,1,2,2,4,7,7,8,10,10,10,10,10]
for i, val in enumerate(knots_maxp):
    knots_maxp[i] = val / 10.0
# control points
ctrl_pts_maxp_2d = [[0,0],
                    [1,1],
                    [2, -1],
                    [3, 0],
                    [5, -2],
                    [6, -3],
                    [8, -4],
                    [10, -5],
                    [10, -6],
                    [9, -6],
                    [8, -5],
                    [7, -4],
                    [6, -2],
                    [5, 0],
                    [3,-1]]

ctrl_pts_maxp_3d = [[0,0,0],
                    [1,1,1],
                    [2, -1, 2],
                    [3, 0, 2],
                    [5, -2, 1],
                    [6, -3, -1],
                    [8, -4, -1],
                    [10, -5, -2],
                    [10, -6, -2],
                    [9, -6, -1],
                    [8, -5, 0],
                    [7, -4, 0],
                    [6, -2, -1],
                    [5, 0, -3],
                    [3,-1, -1]]

ctrl_weights_maxp = [1.0, 2.0, 3.0, 2.0, 1.0,
                             2.0, 3.0, 2.0, 1.0, 2.0, 3.0, 2.0, 1.0, 2.0, 3.0]

# number of knots
m_max = len(knots_maxp)
# number of parameter points
num_points = 100
u = np.linspace(0.0, 1.0, num_points)


def getKnots(p):

    knots = knots_maxp[(pmax - p) : (m_max - (pmax - p))]
    return knots

def getPoints(p, d):

    m = 8 + 2*p
    np = m - p
    # points = [[0 for _ in range(d)] for _ in range(mp)]
    if d == 2:
        points = ctrl_pts_maxp_2d[0:np]
    elif d == 3:
        points = ctrl_pts_maxp_3d[0:np]

    return points

def getWeights(p):

    m = 8 + 2*p
    np = m - p
    weights = ctrl_weights_maxp[0:np]
    return weights



def saveParam(data_out: dict):

    data_out["u"] = dict()
    data_out["u"]["description"] = "This is a set of parameter values at which the various test quantities are computed" 

    data_out["u"]["values"] = u.tolist()


def saveKnots(data_out: dict):

    # ...................................... save knots

    for p in range(1, pmax + 1):

        knots = getKnots(p)
        m = len(knots)
        knot_dataset_str = 'knots_p%i' % p
        data_out[knot_dataset_str] = dict()
        data_out[knot_dataset_str]["description"] = 'Order {} knots'.format(p) 
        data_out[knot_dataset_str]["values"] = knots

def saveWeights(data_out: dict):

    # ...................................... save knots

    for p in range(1, pmax + 1):

        weights = getWeights(p)
        m = len(weights)
        dataset_str = 'weights_p%i' % p
        data_out[dataset_str] = dict()
        data_out[dataset_str]["description"] = "Order {} weights".format(p)
        data_out[dataset_str]["values"] = weights

def saveCtrlpoints(data_out: dict):

    # ...................................... save knots

    for d in range(2, 4):
        for p in range(1, pmax + 1):

            points = np.array(getPoints(p, d))
            dataset_str = 'cpoints_d%i_p%i' % (d, p)
            data_out[dataset_str] = dict()
            data_out[dataset_str]["description"] = "Points for curve of dimension {} and order {}".format(d, p)
            data_out[dataset_str]["values"] = points.tolist()

def curveVals(data_out: dict):

    # .................................. save points
    #set diemnsion
    for d in range(2, 4):
        # Set degree
        for p in range(1, pmax + 1):
            knots = getKnots(p)
            points = getPoints(p, d)
            weights = getWeights(p)
            # create the curve
            curve = NURBS.Curve()
            curve.degree = p
            curve.ctrlpts = points
            curve.knotvector = knots
            curve.weights = weights


            # convert parameters to list
            u_list = u.tolist()
            P = np.array(curve.evaluate_list(u_list))


            dataset_str = 'points_d%i_p%i' % (d, p)
            data_out[dataset_str] = dict()
            data_out[dataset_str]["description"] = "Curve values for dim {} and order {}".format(d, p)
            data_out[dataset_str]["values"] = P.tolist()


def curveDers(data_out: dict):
    #set diemnsion
    for d in range(2, 4):
        # Set degree
        for p in range(1, pmax + 1):
            knots = getKnots(p)
            points = getPoints(p, d)
            weights = getWeights(p)
            # create the curve
            curve = NURBS.Curve()
            curve.degree = p
            curve.ctrlpts = points
            curve.knotvector = knots
            curve.weights = weights

            ders = np.zeros((num_points , (p+1) * d))

            for i, ui in enumerate(u):
                ders_i_tmp = curve.derivatives(ui, p)
                ders_i = list(chain.from_iterable(ders_i_tmp))
                ders[i, :] = ders_i

            dataset_str = 'ders_d%i_p%i' % (d, p)
            data_out[dataset_str] = dict()
            data_out[dataset_str]["description"] = "Der"
            data_out[dataset_str]["values"] = ders.tolist()


def curveTangent(data_out: dict):
    #set diemnsion
    for d in range(2, 4):
        # Set degree
        for p in range(1, pmax + 1):

            knots = getKnots(p)
            points = getPoints(p, d)
            weights = getWeights(p)
            # create the curve
            curve = NURBS.Curve()
            curve.degree = p
            curve.ctrlpts = points
            curve.knotvector = knots
            curve.weights = weights

            tangents = np.zeros((num_points, d))
            tan_i = [0 for _ in range(d)]
            for i, ui in enumerate(u):
                tan_i = operations.tangent(curve, ui, normalize = False)
                tangents[i, :] = tan_i[1]


            dataset_str = 'tangent_d%i_p%i' % (d, p)
            data_out[dataset_str] = dict()
            data_out[dataset_str]["description"] = "tangents of the curve evaluated at the paramters"
            data_out[dataset_str]["values"] = tangents.transpose().tolist()


def curveNormal(data_out: dict):

    #set diemnsion
    for d in range(2, 4):
        # Set degree
        for p in range(1, pmax + 1):

            knots = getKnots(p)
            points = getPoints(p, d)
            weights = getWeights(p)
            # create the curve
            curve = NURBS.Curve()
            curve.degree = p
            curve.ctrlpts = points
            curve.knotvector = knots
            curve.weights = weights

            normals = np.zeros((num_points, 3))
            tan_i = [0 for _ in range(3)]
            for i, ui in enumerate(u):
                tan_i = operations.normal(curve, ui, normalize = False)
                normals[i, :] = tan_i[1]


            dataset_str = 'normal_d%i_p%i' % (d, p)
            data_out[dataset_str]["description"] = "Curve normals evaluated at the parameter values"
            data_out[dataset_str]["values"] = normals.transpose()

# def curveBinormal():

#     f = h5py.File(file_name, 'a')
#     grp = f.create_group('binormals')

#     #set diemnsion
#     for d in range(2, 4):
#         # Set degree
#         for p in range(1, pmax + 1):

#             knots = getKnots(p)
#             points = getPoints(p, d)
#             weights = getWeights(p)
#             # create the curve
#             curve = NURBS.Curve()
#             curve.degree = p
#             curve.ctrlpts = points
#             curve.knotvector = knots
#             curve.weights = weights

#             tangents = np.zeros((num_points, 3))
#             tan_i = [0 for _ in range(3)]
#             for i, ui in enumerate(u):
#                 tan_i = operations.binormal(curve, ui, normalize = False)
#                 tangents[i, :] = tan_i[1]


#             dataset_str = 'binormal_d%i_p%i' % (d, p)
#             dataset = grp.create_dataset(dataset_str,
#                                         (tangents.shape[1], tangents.shape[0]),
#                                         dtype = np.double)
#             dataset[:,:] = tangents.transpose()

def knotInsertion():

    f = h5py.File(file_name, 'a')
    grp = f.create_group('knot_insertion')

    p = 4;
    d = 3

    # ..........................  fisrt insertion
    knots = getKnots(p)
    points = getPoints(p, d)
    weights = getWeights(p)
    # create the curve
    curve1 = NURBS.Curve()
    curve1.degree = p
    curve1.ctrlpts = points
    curve1.knotvector = knots
    curve1.weights = weights

    operations.insert_knot(curve1, [0.5], [4])
    knots1 = np.array(curve1.knotvector)
    points1 = np.array(curve1.ctrlpts)
    pointsw1 = np.array(curve1.ctrlptsw)
    weights1 = np.array(curve1.weights)

    print(knots1)

    grp1 = grp.create_group('insertion1')

    dset = grp1.create_dataset('knot', (1,1), dtype = np.double)
    dset[:,:] = 0.5

    dset = grp1.create_dataset('r', (1,1), dtype = np.int)
    dset[:,:] = 4

    dset = grp1.create_dataset('knots',
                                (1, knots1.shape[0]),
                                dtype = np.double)
    dset[:,:] = knots1.transpose()

    dset = grp1.create_dataset('points', (points1.shape[1], points1.shape[0]),
                                                dtype = np.double)
    dset[:,:] = points1.transpose()

    dset = grp1.create_dataset('pointsw', (pointsw1.shape[1], pointsw1.shape[0]),
                                                dtype = np.double)

    dset[:,:] = pointsw1.transpose()

    dset = grp1.create_dataset('weights', (1, weights1.shape[0]),
                                                dtype = np.double)

    dset[:,:] = weights1.transpose()



    curve2 = NURBS.Curve()
    curve2.degree = p
    curve2.ctrlpts = points
    curve2.knotvector = knots
    curve2.weights = weights

    operations.insert_knot(curve2, [0.7], [2])

    knots2 = np.array(curve2.knotvector)
    points2 = np.array(curve2.ctrlpts)
    pointsw2 = np.array(curve2.ctrlptsw)
    weights2 = np.array(curve2.weights)

    grp1 = grp.create_group('insertion2')

    dset = grp1.create_dataset('knot', (1,1), dtype = np.double)
    dset[:,:] = 0.7

    dset = grp1.create_dataset('r', (1,1), dtype = np.double)
    dset[:,:] = 2

    dset = grp1.create_dataset('knots',
                                (1, knots2.shape[0]),
                                dtype = np.double)
    dset[:,:] = knots2.transpose()

    dset = grp1.create_dataset('points', (points2.shape[1], points2.shape[0]),
                                                dtype = np.double)
    dset[:,:] = points2.transpose()

    dset = grp1.create_dataset('pointsw', (pointsw2.shape[1], pointsw2.shape[0]),
                                                dtype = np.double)

    dset[:,:] = pointsw2.transpose()

    dset = grp1.create_dataset('weights', (1, weights2.shape[0]),
                                                dtype = np.double)

    dset[:,:] = weights2.transpose()



def curveSplit():

    p = 4
    d = 3
    knots = getKnots(p)
    points = getPoints(p, d)
    weights = getWeights(p)

    # -------------------------------------------------------- insertion 1
    # create the curve
    curve = NURBS.Curve()
    curve.degree = p
    curve.ctrlpts = points
    curve.knotvector = knots
    curve.weights = weights


    u = 0.5
    curves = operations.split_curve(curve, u)

    knots1 = np.array(curves[0].knotvector)
    points1 = np.array(curves[0].ctrlpts)
    pointsw1 = np.array(curves[0].ctrlptsw)
    weights1 = np.array(curves[0].weights)

    knots2 = np.array(curves[1].knotvector)
    points2 = np.array(curves[1].ctrlpts)
    pointsw2 = np.array(curves[1].ctrlptsw)
    weights2 = np.array(curves[1].weights)



    f = h5py.File(file_name, 'a')
    grp = f.create_group('curve_split')

    grp0 = grp.create_group('split1/curve0')
    grp1 = grp.create_group('split1/curve1')

    dset = grp0.create_dataset('knots',
                                (1, knots1.shape[0]),
                                dtype = np.double)
    dset[:,:] = knots1.transpose()


    dset = grp0.create_dataset('points',
                                (points1.shape[1], points1.shape[0]),
                                dtype = np.double)
    dset[:,:] = points1.transpose()

    dset = grp0.create_dataset('pointsw',
                                (pointsw1.shape[1], pointsw1.shape[0]),
                                dtype = np.double)
    dset[:,:] = pointsw1.transpose()

    dset = grp0.create_dataset('weights',
                                (1, weights1.shape[0]),
                                dtype = np.double)
    dset[:,:] = weights1.transpose()


    # -------
    dset = grp1.create_dataset('knots',
                                (1, knots2.shape[0]),
                                dtype = np.double)
    dset[:,:] = knots2.transpose()


    dset = grp1.create_dataset('points',
                                (points2.shape[1], points2.shape[0]),
                                dtype = np.double)
    dset[:,:] = points2.transpose()

    dset = grp1.create_dataset('pointsw',
                                (pointsw2.shape[1], pointsw2.shape[0]),
                                dtype = np.double)
    dset[:,:] = pointsw2.transpose()

    dset = grp1.create_dataset('weights',
                                (1, weights2.shape[0]),
                                dtype = np.double)
    dset[:,:] = weights2.transpose()

    # -------------------------------------------------------- insertion 2
    # create the curve
    curve = NURBS.Curve()
    curve.degree = p
    curve.ctrlpts = points
    curve.knotvector = knots
    curve.weights = weights


    u = 0.7
    curves = operations.split_curve(curve, u)

    knots1 = np.array(curves[0].knotvector)
    points1 = np.array(curves[0].ctrlpts)
    pointsw1 = np.array(curves[0].ctrlptsw)
    weights1 = np.array(curves[0].weights)

    knots2 = np.array(curves[1].knotvector)
    points2 = np.array(curves[1].ctrlpts)
    pointsw2 = np.array(curves[1].ctrlptsw)
    weights2 = np.array(curves[1].weights)


    grp = f['curve_split']
    grp0 = grp.create_group('split2/curve0')
    grp1 = grp.create_group('split2/curve1')

    dset = grp0.create_dataset('knots',
                                (1, knots1.shape[0]),
                                dtype = np.double)
    dset[:,:] = knots1.transpose()


    dset = grp0.create_dataset('points',
                                (points1.shape[1], points1.shape[0]),
                                dtype = np.double)
    dset[:,:] = points1.transpose()

    dset = grp0.create_dataset('pointsw',
                                (pointsw1.shape[1], pointsw1.shape[0]),
                                dtype = np.double)
    dset[:,:] = pointsw1.transpose()

    dset = grp0.create_dataset('weights',
                                (1, weights1.shape[0]),
                                dtype = np.double)
    dset[:,:] = weights1.transpose()


    # -------
    dset = grp1.create_dataset('knots',
                                (1, knots2.shape[0]),
                                dtype = np.double)
    dset[:,:] = knots2.transpose()


    dset = grp1.create_dataset('points',
                                (points2.shape[1], points2.shape[0]),
                                dtype = np.double)
    dset[:,:] = points2.transpose()

    dset = grp1.create_dataset('pointsw',
                                (pointsw2.shape[1], pointsw2.shape[0]),
                                dtype = np.double)
    dset[:,:] = pointsw2.transpose()

    dset = grp1.create_dataset('weights',
                                (1, weights2.shape[0]),
                                dtype = np.double)
    dset[:,:] = weights2.transpose()



def saveIntegration():

    p = 3
    knots = [0,0,0,0, 0.22, 0.55, 0.55, 0.8, 1,1,1,1]
    points = [[0,0,0],
              [1,0,0],
              [2,1,1],
              [3,4,2],
              [4,5,3],
              [5,3,2],
              [6,4,1],
              [7,4,0]]


    weights = [1,1,1,2,1.5,1,1,0.7]

    curve = NURBS.Curve()
    curve.degree = p
    curve.ctrlpts = points
    curve.knotvector = knots
    curve.weights = weights


    knots1 = np.array(curve.knotvector)
    points1 = np.array(curve.ctrlpts)
    pointsw1 = np.array(curve.ctrlptsw)
    weights1 = np.array(curve.weights)


    def jacobian(u):
        tan_i = operations.tangent(curve, u, normalize = False)
        jac = math.sqrt(np.dot(tan_i[1], tan_i[1]))
        return jac

    def fcn3d(x, y, z):
        val = x * y * z * (math.exp(y)/500)
        return val

    def integrand(u):
        X = curve.evaluate_single(u)
        x = X[0]
        y = X[1]
        z = X[2]
        jac = jacobian(u)
        f = jac * fcn3d(x,y,z)
        return f


    If, err = sci.quad(integrand, 0, 1, epsabs = 1e-6, epsrel = 1e-6)

    f = h5py.File(file_name, 'a')
    grp = f.create_group('integration')

    dset = grp.create_dataset('p', (1,0), dtype = np.int)
    dset[0] = p


    dset = grp.create_dataset('knots',
                                (1, knots1.shape[0]),
                                dtype = np.double)
    dset[:,:] = knots1.transpose()

    dset = grp.create_dataset('points', (points1.shape[1], points1.shape[0]),
                                                dtype = np.double)
    dset[:,:] = points1.transpose()

    dset = grp.create_dataset('pointsw', (pointsw1.shape[1], pointsw1.shape[0]),
                                                dtype = np.double)

    dset[:,:] = pointsw1.transpose()

    dset = grp.create_dataset('weights', (1, weights1.shape[0]),
                                                dtype = np.double)

    dset[:,:] = weights1.transpose()

    dset = grp.create_dataset('integral', (1,), dtype = np.double)
    dset[0] = If

    # curve.vis = VisMPL.VisCurve3D()
    # curve.render()

    #
    #
    # n = 1000
    # delta_u = 1 / n
    # u = [delta_u*i for i in range(0,n+1)]
    # jac = [0 for i in range(0,n+1)]
    # fcn_vals = [0 for i in range(0,n+1)]
    #
    # for i, ui in enumerate(u):
    #     jac[i] = jacobian(ui)
    #     fcn_vals[i] = integrand(ui)
    #
    # #
    # fig, ax = plt.subplots()
    # ax.plot(u, fcn_vals)
    # plt.show()








def main():

    data_out = dict()
    saveParam(data_out)
    saveKnots(data_out)
    saveWeights(data_out)
    saveCtrlpoints(data_out)
    curveVals(data_out)
    curveDers(data_out)
    curveTangent(data_out)

    with open(file_name, 'w') as f:
        json.dump(data_out, f, indent=4)
    # curveNormal(data_out)
    # curveBinormal(data_out)
    # knotInsertion(data_out)
    # curveSplit(data_out)
    # saveIntegration(data_out)




if __name__ == '__main__':
    main()
