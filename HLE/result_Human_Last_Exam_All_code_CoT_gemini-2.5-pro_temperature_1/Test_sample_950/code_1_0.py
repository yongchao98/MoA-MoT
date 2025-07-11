def solve_grassmannian_torsion_rank():
    """
    Calculates the rank of the torsion subgroup of the integral cohomology
    ring of the space of 3-subspaces of R^5, i.e., Gr(3, 5).
    """
    k = 3
    n = 5
    n_minus_k = n - k

    partitions = []
    # Generate all partitions lambda = (l1, l2, l3) such that
    # n-k >= l1 >= l2 >= l3 >= 0.
    for l1 in range(n_minus_k, -1, -1):
        for l2 in range(l1, -1, -1):
            for l3 in range(l2, -1, -1):
                partitions.append((l1, l2, l3))

    # Total number of Schubert cells, which is the sum of Z_2-Betti numbers.
    sum_d_i = len(partitions)

    # Count partitions where all parts are even. This gives the sum of Betti numbers.
    even_partitions = []
    for p in partitions:
        if all(val % 2 == 0 for val in p):
            even_partitions.append(p)
    
    sum_beta_i = len(even_partitions)

    # The total rank of the torsion subgroup is (sum_d_i - sum_beta_i) / 2
    torsion_rank = (sum_d_i - sum_beta_i) // 2

    print(f"The space is the real Grassmannian Gr({k}, {n}).")
    print(f"The Schubert cells are indexed by partitions lambda = (l1, ..., l{k}) fitting in a {k}x{n_minus_k} box.")
    print(f"List of all partitions: {partitions}")
    print(f"The sum of Z2-Betti numbers is the total number of partitions.")
    print(f"Sum(d_i) = {sum_d_i}")
    print("-" * 20)
    print("The sum of integral Betti numbers is the number of partitions where all parts are even.")
    print(f"List of partitions with even parts: {even_partitions}")
    print(f"Sum(beta_i) = {sum_beta_i}")
    print("-" * 20)
    print("The total rank of the torsion subgroup T is given by the formula:")
    print("T = (Sum(d_i) - Sum(beta_i)) / 2")
    print(f"T = ({sum_d_i} - {sum_beta_i}) / 2")
    print(f"T = {sum_d_i - sum_beta_i} / 2")
    print(f"T = {torsion_rank}")
    print("-" * 20)
    print("The rank of the torsion subgroup of the integral cohomology ring is:")
    print(torsion_rank)

solve_grassmannian_torsion_rank()