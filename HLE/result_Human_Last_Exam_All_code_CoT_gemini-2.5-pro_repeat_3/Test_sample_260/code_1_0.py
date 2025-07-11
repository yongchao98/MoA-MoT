from fractions import Fraction

def solve_rank_of_kernel():
    """
    Calculates the rank of the kernel of the Hurewicz homomorphism for Y = X1 # X2 # X3.
    """
    # Step 1: Define the orders of the fundamental groups of X1, X2, X3.
    # pi_1(X1) = Z_5, pi_1(X2) = Z_8, pi_1(X3) = Z_2.
    d1 = 5
    d2 = 8
    d3 = 2
    orders = [d1, d2, d3]
    n = len(orders)

    print("Step-by-step calculation of the rank 'r'.")
    print("-" * 40)
    print(f"The fundamental groups of the spaces are Z_{d1}, Z_{d2}, and Z_{d3}.")
    print(f"The fundamental group of the connected sum is G = Z_{d1} * Z_{d2} * Z_{d3}.")
    print(f"The first homology group is H1(G) = Z_{d1} + Z_{d2} + Z_{d3}.")
    print("The kernel K of the Hurewicz map is the commutator subgroup [G,G], which is a free group.")
    print("The rank 'r' of K is computed using the formula: r = 1 - |H1(G)| * ( sum(1/|Gi|) - (n-1) )")
    print("-" * 40)

    # Step 2: Calculate the order of the first homology group H1(Y).
    # H1(Y) = Z_5 + Z_8 + Z_2.
    # |H1(Y)| = 5 * 8 * 2 = 80.
    h1_order = d1 * d2 * d3

    # Step 3: Calculate the Euler characteristic of the group G = pi_1(Y).
    # chi(G) = chi(Z_5) + chi(Z_8) + chi(Z_2) - (3-1)
    # chi(Z_d) = 1/d
    chi_G_sum_inv = sum(Fraction(1, d) for d in orders)
    chi_G = chi_G_sum_inv - (n - 1)

    # Step 4: Calculate the Euler characteristic of the kernel K.
    # chi(K) = chi(G) * |H1(Y)|
    chi_K = chi_G * h1_order

    # Step 5: Calculate the rank r from chi(K) = 1 - r.
    rank = 1 - chi_K
    
    # Print the detailed equation
    sum_inv_str = " + ".join([f"1/{d}" for d in orders])
    prod_orders_str = " * ".join([str(d) for d in orders])
    print("The final equation for the rank 'r' is:")
    print(f"r = 1 - ( ({sum_inv_str}) - ({n}-1) ) * ({prod_orders_str})")
    print("\nCalculating the components of the equation:")
    print(f"1. Sum of reciprocal orders: {sum_inv_str} = {float(chi_G_sum_inv):.3f}")
    print(f"2. Euler characteristic of G: chi(G) = {float(chi_G_sum_inv):.3f} - {n-1} = {float(chi_G)}")
    print(f"3. Order of H1(G): |H1(G)| = {prod_orders_str} = {h1_order}")
    print(f"4. Euler characteristic of K: chi(K) = chi(G) * |H1(G)| = {float(chi_G)} * {h1_order} = {int(chi_K)}")
    print(f"5. Rank of K: r = 1 - chi(K) = 1 - ({int(chi_K)}) = {int(rank)}")

    print("\nFinal Answer:")
    print(f"The rank of the kernel K is {int(rank)}.")

solve_rank_of_kernel()