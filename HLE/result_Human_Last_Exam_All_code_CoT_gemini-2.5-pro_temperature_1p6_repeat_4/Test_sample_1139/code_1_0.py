def solve_class_d_variables():
    """
    Calculates the number of non-Grassmann variables for the supersymmetric sigma-model
    for disordered systems of symmetry class D with two replicas.
    """
    # Number of replicas for the model
    n_r = 2

    print(f"Starting calculation for a class D system with {n_r} replicas.")
    print("-" * 50)

    # Define a function to calculate the dimension of the special orthogonal group SO(m)
    def dim_so(m):
        """Calculates the dimension of the Lie group SO(m)."""
        return m * (m - 1) // 2

    # Define a function to calculate the dimension of the symplectic group Sp(2k)
    def dim_sp(two_k):
        """Calculates the dimension of the Lie group Sp(2k)."""
        k = two_k // 2
        return k * (2 * k + 1)

    # --- Step 1: Calculate the dimension of the numerator group G_0 ---
    print("Calculating dimension of the numerator's bosonic subgroup G_0 = SO(2*n_r) x Sp(2*n_r):")
    g_so_size = 2 * n_r
    g_sp_size = 2 * n_r

    dim_g_so = dim_so(g_so_size)
    dim_g_sp = dim_sp(g_sp_size)
    dim_g_0 = dim_g_so + dim_g_sp
    
    print(f"  dim(SO({g_so_size})) = {g_so_size}*({g_so_size}-1)/2 = {dim_g_so}")
    print(f"  dim(Sp({g_sp_size})) = ({g_sp_size}//2)*(2*({g_sp_size}//2)+1) = {dim_g_sp}")
    print(f"  Dimension of G_0 = {dim_g_so} + {dim_g_sp} = {dim_g_0}")
    print("-" * 50)

    # --- Step 2: Calculate the dimension of the denominator group H_0 ---
    print("Calculating dimension of the denominator's bosonic subgroup H_0 = (SO(n_r) x Sp(n_r))^2:")
    h_so_size = n_r
    h_sp_size = n_r

    dim_h_so_factor = dim_so(h_so_size)
    dim_h_sp_factor = dim_sp(h_sp_size)
    dim_h_0_factor = dim_h_so_factor + dim_h_sp_factor
    dim_h_0 = 2 * dim_h_0_factor

    print(f"  Dimension of one factor SO({h_so_size}) x Sp({h_sp_size}):")
    print(f"    dim(SO({h_so_size})) = {h_so_size}*({h_so_size}-1)/2 = {dim_h_so_factor}")
    print(f"    dim(Sp({h_sp_size})) = ({h_sp_size}//2)*(2*({h_sp_size}//2)+1) = {dim_h_sp_factor}")
    print(f"    Subtotal = {dim_h_so_factor} + {dim_h_sp_factor} = {dim_h_0_factor}")
    print(f"  Total dimension of H_0 = 2 * {dim_h_0_factor} = {dim_h_0}")
    print("-" * 50)

    # --- Step 3: Final Calculation ---
    print("The number of variables is dim(G_0) - dim(H_0).")
    num_variables = dim_g_0 - dim_h_0
    print(f"Final equation: {dim_g_0} - {dim_h_0} = {num_variables}")
    print("-" * 50)
    print(f"\nThe required number of non-Grassmann variables is {num_variables}.")

solve_class_d_variables()