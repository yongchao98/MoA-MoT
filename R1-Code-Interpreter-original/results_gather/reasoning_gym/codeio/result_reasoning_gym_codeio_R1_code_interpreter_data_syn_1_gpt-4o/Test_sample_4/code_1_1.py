import numpy as np
import math
from itertools import combinations_with_replacement
from itertools import permutations

def how_many_derivs(k, n):
    val = 1
    for i in range(n):
        val *= (3 * k + i)
    return int((1 / math.factorial(n)) * val)

def generate_buffer_lookup(nparams, deriv_order):
    dimensions = (nparams,) * deriv_order 
    buffer_index_lookup = np.zeros(dimensions, dtype=int)
    count = 0
    for idx in combinations_with_replacement(np.arange(nparams), deriv_order):
        for perm in permutations(idx):
            buffer_index_lookup[perm] = count
        count += 1
    return buffer_index_lookup

def generate_multi_index_lookup(nparams, deriv_order, nderivs):
    lookup = np.zeros((nderivs, deriv_order), int)
    idx = 0
    for indices in combinations_with_replacement(np.arange(nparams), deriv_order):
        for i in range(len(indices)):
            lookup[idx, i] = indices[i]
        idx += 1
    return lookup

def main_solution(ncenters, deriv_order):
    ncenters = int(ncenters)
    deriv_order = int(deriv_order)
    
    nderivs = how_many_derivs(ncenters, deriv_order)
    nparams = ncenters * 3
    
    lookup_forward = generate_multi_index_lookup(nparams, deriv_order, nderivs)
    lookup_backward = generate_buffer_lookup(nparams, deriv_order)
    
    return {
        "lookup_forward_shape": lookup_forward.shape,
        "lookup_backward_shape": lookup_backward.shape
    }

# Test the function with the inferred inputs
result = main_solution(1, 1)
print(result)