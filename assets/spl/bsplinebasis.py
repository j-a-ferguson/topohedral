
from geomdl import helpers
import numpy as np
import json, os

file_name = 'bsplinebasis.json'
# max order
pmax = 4
# Set knot vector
knots_maxp = [0, 0, 0, 0, 0, 1, 2, 4, 4,  5, 6, 6, 6, 6, 6]


# max len
m_max = len(knots_maxp)
# number of parameter points
num_points = 100
u = np.linspace(0.0, 6.0, num_points)

def getKnots(p):

    knots = knots_maxp[(pmax - p) : (m_max - (pmax - p))]
    return knots


def saveParam(data_out: dict):

    data_out["u"] = dict()
    data_out["u"]["description"] = '''100 Parameter values in range [0, 6]'''    
    data_out["u"]["values"] = u.tolist()



def saveKnots(data_out: dict):

    # ...................................... save knots

    for p in range(pmax + 1):

        knots = getKnots(p)
        m = len(knots)
        knot_dataset_str = 'knots_p%i' % p
        data_out[knot_dataset_str] = dict()
        data_out[knot_dataset_str]["descrption"] = 'knots for order %d' % (p)
        data_out[knot_dataset_str]["values"] = knots


def spans(data_out: dict):

    spans = np.zeros((num_points), dtype = np.uint64)
    for p in range(pmax + 1):
        knots_u = getKnots(p)
        m = len(knots_u)
        num_basis = m - p - 1

        for i, ui in enumerate(u):

            spans[i] = helpers.find_span_binsearch(p, knots_u, num_basis, ui)

        span_dataset_str = 'span_p%i' % (p)
        data_out[span_dataset_str] = dict()
        data_out[span_dataset_str]["description"] = 'knot span value at every parameter value'
        data_out[span_dataset_str]["values"] = spans.tolist()


def basisFuns(data_out: dict):
    for p in range(0,pmax + 1):

        knots_u = getKnots(p)
        m = len(knots_u)
        num_basis = m - p - 1

        Njp = np.zeros((num_points, num_basis))

        for i, ui in enumerate(u):

            span = helpers.find_span_binsearch(p, knots_u, num_basis, ui)
            vals = helpers.basis_function(p, knots_u, span, ui)
            # print('%i %f %i' % (i, ui, span))
            for j in range(span - p, span + 1):
                Njp[i, j] = vals[j - (span - p)]


        #........................................ save basis
        basis_dataset_str = 'basis_p%i' % p
        data_out[basis_dataset_str] = dict()
        data_out[basis_dataset_str]["description"] = 'shape functions for order %d' % (p)
        data_out[basis_dataset_str]["values"] =  Njp.tolist()



def basisDers(data_out: dict):
    for p in range(pmax+1):
        k = p
        knots_u = getKnots(p)
        m = len(knots_u)
        num_basis = m - p - 1

        dNjp = np.zeros((num_points * (k+1), num_basis))

        for i, ui in enumerate(u):
            start = i * (k+1)
            end = (i + 1) * (k+1)

            span = helpers.find_span_binsearch(p, knots_u, num_basis, ui)
            vals_i = np.array(helpers.basis_function_ders(p, knots_u, span, ui, k))

            dNjp[start:end, (span-p):(span+1)] = vals_i


        dataset_str = 'ders_p%i' % p
        data_out[dataset_str] = dict()
        data_out[dataset_str]["description"] = 'shape function derivatives for order %d' % (p)
        data_out[dataset_str]["values"] = dNjp.tolist()

def main():


    data_out = dict()
    saveParam(data_out)
    saveKnots(data_out)
    spans(data_out)
    basisFuns(data_out)
    basisDers(data_out)

    print(data_out)

    with open(file_name, 'w') as f:
        json.dump(data_out, f, indent=4)








if __name__ == '__main__':
    main()
