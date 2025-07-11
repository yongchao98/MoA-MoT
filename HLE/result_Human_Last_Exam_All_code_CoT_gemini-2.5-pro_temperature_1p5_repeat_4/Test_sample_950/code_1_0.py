import numpy as np

def solve_cohomology_torsion_rank():
    """
    Calculates the rank of the torsion subgroup of the integral cohomology ring
    of the space of 3-subspaces of R^5.
    """
    # Step 1: Identify the space
    # The space of 3-subspaces of R^5 is the Grassmannian G(3,5).
    # G(k,n) is homeomorphic to G(n-k, n). So G(3,5) is same as G(2,5).
    k = 2
    n = 5
    dim = k * (n - k)

    print(f"The space is the real Grassmannian G_3(R^5), which is homeomorphic to G_{k}(R^{n}) for k={k}, n={n}.")
    print(f"The dimension of this manifold is k*(n-k) = {dim}.")
    print("-" * 30)

    # Step 2: Betti numbers with Z_2 coefficients
    # The Z_2-Poincare polynomial is the q-binomial coefficient (n k)_t.
    # For G(2,5), this is (1-t^5)(1-t^4)/((1-t)(1-t^2)) = (1+t+t^2+t^3+t^4)(1+t^2)
    p1 = np.poly1d([1, 1, 1, 1, 1])
    p2 = np.poly1d([1, 0, 1])
    p_z2 = p1 * p2
    d_coeffs = p_z2.coeffs.astype(int)
    
    # Pad with zeros if necessary for all dimensions
    d = np.zeros(dim + 1, dtype=int)
    d[:len(d_coeffs)] = d_coeffs
    
    print("Step 2: Calculate Z_2-Betti numbers (d_i)")
    print(f"The Z_2-Poincare polynomial P(t) = {p_z2}")
    print(f"The coefficients d_i = dim(H^i(G; Z_2)) are:")
    for i, val in enumerate(d):
        print(f"d_{i} = {val}")
    print("-" * 30)

    # Step 3: Betti numbers with Q coefficients
    # The rational Poincare polynomial for G(2,5) is (1-t^8)/(1-t^4) = 1 + t^4
    b = np.zeros(dim + 1, dtype=int)
    b[0] = 1
    b[4] = 1

    print("Step 3: Calculate rational Betti numbers (b_i)")
    print("The rational Poincare polynomial is P(t) = 1 + t^4.")
    print("The coefficients b_i = rank(H^i(G; Z)) are:")
    for i, val in enumerate(b):
        print(f"b_{i} = {val}")
    print("-" * 30)

    # Step 4: Determine the torsion ranks t_i
    # Using the formula from UCT: d_i = b_i + t_i + t_{i+1}
    # This implies t_i = d_i - b_i - t_{i+1}
    t = np.zeros(dim + 2, dtype=int)  # t_0, ..., t_dim, t_{dim+1}

    print("Step 4: Calculate torsion subgroup ranks (t_i)")
    print("Using the recursive formula t_i = d_i - b_i - t_{i+1}, starting from t_{dim+1}=0.")
    for i in range(dim, -1, -1):
        t[i] = d[i] - b[i] - t[i+1]
        print(f"t_{i} = d_{i} - b_{i} - t_{i+1} = {d[i]} - {b[i]} - {t[i+1]} = {t[i]}")
    
    t_list = t[0:dim+1].tolist()
    print("-" * 30)

    # Step 5: Calculate the total rank
    total_rank = np.sum(t)
    
    print("Step 5: Calculate the total rank of the torsion subgroup")
    print(f"The ranks of the torsion subgroups are: {t_list}")
    
    non_zero_ti = {i: val for i, val in enumerate(t_list) if val > 0}
    sum_str = " + ".join([f"t_{i}" for i in non_zero_ti.keys()])
    val_str = " + ".join([str(val) for val in non_zero_ti.values()])
    
    print(f"The total rank is the sum: {sum_str} = {val_str} = {total_rank}")
    print("-" * 30)

    return total_rank

if __name__ == '__main__':
    final_answer = solve_cohomology_torsion_rank()
    # The format required by the calling system.
    # print(f"\nFinal Answer: {final_answer}")
    
solve_cohomology_torsion_rank()