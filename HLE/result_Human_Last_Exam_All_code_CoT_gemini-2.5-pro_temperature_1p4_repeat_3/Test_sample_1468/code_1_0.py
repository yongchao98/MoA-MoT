import sympy

def solve_lower_bound():
    """
    This function demonstrates the dimensional argument for the lower bound of m.
    """
    # Define symbolic variables for the dimensions involved.
    # N: number of data points
    # d_prime: dimension of the z vectors
    # q: sparsity
    # m: width of the hidden layer (the variable we want to bound)
    N, d_prime, q, m = sympy.symbols('N d_prime q m', integer=True, positive=True)

    # The proof focuses on the network's ability to process a single z_k vector.
    # The components for z_k live in a d'-dimensional subspace of the full input space.
    # Let's call this subspace V_k_z.
    dim_V_k_z = d_prime

    # The total dimension of the flattened input x is N * d, where d = d' + q + 1.
    d = d_prime + q + 1
    dim_total = N * d

    # The model's first layer is a linear transformation W of size m x (Nd).
    # The kernel (or null space) of W is the set of input vectors that W maps to zero.
    # By the rank-nullity theorem, dim(ker(W)) = dim_total - rank(W).
    # Since W has m rows, its rank can be at most m. So, rank(W) <= m.
    # This gives a lower bound on the dimension of the kernel.
    dim_ker_W_lower_bound = dim_total - m

    # The core of the argument is to analyze the intersection of the subspace V_k_z
    # and the kernel of W, ker(W). The dimension of the intersection of two subspaces
    # U and V is bounded by: dim(U ∩ V) >= dim(U) + dim(V) - dim(Total Space).
    dim_intersection_lower_bound = dim_V_k_z + dim_ker_W_lower_bound - dim_total

    # Let's substitute our expressions into this formula.
    dim_intersection_formula = d_prime + (N * d - m) - (N * d)

    # Sympy can simplify this for us.
    final_dim_intersection = sympy.simplify(dim_intersection_formula)

    print("### Proof via Dimensional Analysis ###")
    print(f"Let's analyze the input space and the effect of the weight matrix W.")
    print(f"1. The input vector z_k is {d_prime}-dimensional. The corresponding components in the flattened input x span a subspace V_k_z of dimension d' = {d_prime}.")
    print(f"2. The weight matrix W is size {m} x {N*d}. Its kernel, ker(W), has dimension at least Nd - m = {dim_ker_W_lower_bound}.")
    print(f"3. We consider the intersection of these two subspaces, V_k_z and ker(W).")
    print(f"   dim(V_k_z ∩ ker(W)) >= dim(V_k_z) + dim(ker(W)) - dim(Total Space)")
    print(f"   dim(V_k_z ∩ ker(W)) >= {d_prime} + ({N*d - m}) - {N*d}")
    print(f"   Simplifying this gives: dim(V_k_z ∩ ker(W)) >= {final_dim_intersection}")
    
    print("\n4. If m < d', then d' - m > 0. This means the intersection is non-trivial.")
    print("   This implies there exists a non-zero vector Δz (a change in z_k) such that W(Δz) = 0.")
    print("   The network is 'blind' to this change.")
    print("\n5. If the network is blind to Δz, it produces the same output for an input x and x+Δz.")
    print("   However, the target qSA value changes. For appropriately chosen inputs, ||qSA(x) - qSA(x+Δz)|| can be as large as 2/q.")
    print("   The network error condition requires ||qSA(x) - qSA(x+Δz)|| <= 2ε, which is 1/q.")
    print("   This leads to the contradiction 2/q <= 1/q, which is false.")

    print("\n6. To avoid this contradiction, the intersection must be trivial (dimension 0).")
    print("   This requires the lower bound on its dimension to be less than or equal to 0.")
    print(f"   So, we must have: {final_dim_intersection} <= 0")
    print(f"   Which leads to the final inequality:")
    print(f"   m >= {d_prime}")


solve_lower_bound()