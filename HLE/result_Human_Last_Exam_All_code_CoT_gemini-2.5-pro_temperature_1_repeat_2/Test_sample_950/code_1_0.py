import numpy as np

def solve():
    """
    Calculates the rank of the torsion subgroup of the integral cohomology ring
    of the space of 3-subspaces of R^5, which is the Grassmannian Gr(3,5).
    """
    # Dimension of the manifold
    k = 3
    n = 5
    dim = k * (n - k)

    # h_i = dim(H^i(Gr(3,5); Z/2Z)) from Schubert calculus
    # Partitions (l_1, l_2, l_3) with 2 >= l_1 >= l_2 >= l_3 >= 0
    # deg 0: (0,0,0) -> 1
    # deg 1: (1,0,0) -> 1
    # deg 2: (2,0,0), (1,1,0) -> 2
    # deg 3: (2,1,0), (1,1,1) -> 2
    # deg 4: (2,2,0), (2,1,1) -> 2
    # deg 5: (2,2,1) -> 1
    # deg 6: (2,2,2) -> 1
    h = np.array([1, 1, 2, 2, 2, 1, 1])
    print(f"The dimensions of the mod 2 cohomology groups h_i are: {h.tolist()}")

    # beta_i = rank(H^i(Gr(3,5); Z)), from rational cohomology calculation.
    # H*(Gr(3,5); Q) = Q. So beta_0=1, beta_i=0 for i>0.
    beta = np.zeros(dim + 1, dtype=int)
    beta[0] = 1
    print(f"The integral Betti numbers beta_i are: {beta.tolist()}")

    # From UCT, h_i = beta_i + t_{i+1}, where t_i = rank(Tors(H^i))
    # We can compute t_i for i=1,...,dim. t_0 is always 0.
    t = np.zeros(dim + 1, dtype=int)
    print("\nUsing the Universal Coefficient Theorem formula h_i = beta_i + t_{i+1}:")
    for i in range(dim):
        t[i+1] = h[i] - beta[i]
        print(f"For i={i}: h_{i} = beta_{i} + t_{i+1}")
        print(f"         {h[i]} = {beta[i]} + t_{i+1}  => t_{i+1} = {t[i+1]}")

    print(f"\nThe ranks of the torsion subgroups T_i are t_i: {t.tolist()}")

    # The rank of the torsion subgroup of the integral cohomology ring
    # is the sum of the ranks of the torsion parts of each group.
    total_torsion_rank = np.sum(t)
    
    equation_parts = [f"t_{i}" for i in range(1, dim + 1)]
    equation_values = [str(val) for val in t[1:]]
    
    print("\nThe total rank of the torsion subgroup is the sum:")
    print(f"Rank = {' + '.join(equation_parts)} = {' + '.join(equation_values)} = {total_torsion_rank}")
    
    print("\nFinal Answer:")
    print(total_torsion_rank)

solve()