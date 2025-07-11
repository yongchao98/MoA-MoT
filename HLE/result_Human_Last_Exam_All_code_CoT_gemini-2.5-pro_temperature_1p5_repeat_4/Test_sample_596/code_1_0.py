import sys

# Plan
# 1. Define the Betti numbers of the Mobius band (K) over the field F_2.
#    The Mobius band is a 2-dimensional simplicial complex.
#    b_0(K) = 1 (it is path-connected).
#    b_1(K) = 1 (it deformation retracts onto a circle).
#    b_2(K) = 0 (it is a manifold with boundary).
betti_K = {0: 1, 1: 1, 2: 0}

# 2. Determine the relative Betti numbers of the pair (K, partial K) where
#    partial K is the boundary of the Mobius band.
#    Using Lefschetz duality for a non-orientable n-manifold over F_2, we have:
#    b_i(K, partial K) = b_{n-i}(K). Here n=2.
#    b_0(K, partial K) = b_2(K) = 0.
#    b_1(K, partial K) = b_1(K) = 1.
#    b_2(K, partial K) = b_0(K) = 1.
betti_K_partialK = {
    0: betti_K[2],
    1: betti_K[1],
    2: betti_K[0]
}

# 3. The number of higher dimensional rooted forests that fail to collapse
#    is given by the formula: Sum_{i=0 to n} b_i(K) * b_i(K, partial K).
term0 = betti_K[0] * betti_K_partialK[0]
term1 = betti_K[1] * betti_K_partialK[1]
term2 = betti_K[2] * betti_K_partialK[2]
total_sum = term0 + term1 + term2

# 4. Print the final calculation, showing each number in the final equation.
# The calculation is sum_{i=0 to 2} b_i(K) * b_i(K, partial K).
# = b_0(K)*b_0(K,pK) + b_1(K)*b_1(K,pK) + b_2(K)*b_2(K,pK)
# = b_0(K)*b_2(K) + b_1(K)*b_1(K) + b_2(K)*b_0(K)
print(f"{betti_K[0]} * {betti_K_partialK[0]} + {betti_K[1]} * {betti_K_partialK[1]} + {betti_K[2]} * {betti_K_partialK[2]} = {total_sum}")