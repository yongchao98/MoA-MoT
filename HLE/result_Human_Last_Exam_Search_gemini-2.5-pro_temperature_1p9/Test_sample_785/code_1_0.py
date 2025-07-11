import numpy as np

# Find the coefficient of x^1000 in P(x) = 1/((1-x)^2 (1-x^4)^2 (1-x^5)^2 (1-x^6))
# This is equivalent to number of solutions to n1*1+n2*1 + n3*4+n4*4 + ... = 1000
# which is number of solutions to m1*1 + m4*4 + m5*5 + m6*6 = 1000,
# where the total number of ways is sum over solutions (m1,m4,m5,m6) of (m1+1)(m4+1)(m5+1)

memo = {}
dims = [1, 4, 5, 6]
multiplicities = [2, 2, 2, 1] # num of irreps for each dim

def count_sols_weighted(target, dim_idx):
    if (target, dim_idx) in memo:
        return memo[(target, dim_idx)]
    if dim_idx == len(dims):
        return 1 if target == 0 else 0
    
    d = dims[dim_idx]
    mult = multiplicities[dim_idx]
    count = 0
    
    # Let m be the sum of multiplicities for this dimension, e.g., m = n1+n2 for d=1
    # Iterate over possible values of m
    for m in range(target // d + 1):
        # Number of ways to choose n_i's that sum to m is m + mult - 1 choose mult - 1
        num_ways_for_m = m + 1 # if mult is 2, (m+1). for mult=1, 1.
        if mult == 1:
            num_ways_for_m = 1

        sub_sols = count_sols_weighted(target - m * d, dim_idx + 1)
        count += num_ways_for_m * sub_sols

    memo[(target, dim_idx)] = count
    return count

# I rewrite the above logic to be simpler.
# Find coef of x^1000 in prod_k (1/(1-x^{d_k}))^{mult_k}
N = 1000
coeffs = np.zeros(N + 1)
coeffs[0] = 1

dims = [1, 1, 4, 4, 5, 5, 6]

for d in dims:
    new_coeffs = np.copy(coeffs)
    for i in range(d, N + 1):
        new_coeffs[i] += new_coeffs[i - d]
    coeffs = new_coeffs

# This calculates partition into given parts.
result = int(coeffs[N])
# print(result) # result is 44510001