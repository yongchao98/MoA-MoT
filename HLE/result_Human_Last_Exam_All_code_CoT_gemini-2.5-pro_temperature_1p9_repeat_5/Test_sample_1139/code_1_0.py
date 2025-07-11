def solve_class_d_variables():
    """
    Calculates the number of non-Grassman variables for the supersymmetric sigma-model
    of class D with a specified number of replicas.
    """
    N_replicas = 2
    
    print("This script calculates the number of non-Grassman (bosonic) variables required to parametrize")
    print("the sigma-model for disordered systems of symmetry class D with 2 replicas.")
    print("-" * 70)
    print(f"Problem parameters: Symmetry Class = D, Number of Replicas (N) = {N_replicas}\n")

    # --- Theory ---
    print("Step 1: Identify the structure of the bosonic manifold.")
    print("For class D, the bosonic sector of the sigma-model is the manifold:")
    print("  M_b = (SO(4N) / (SO(2N) x SO(2N))) x (Sp(4N) / (Sp(2N) x Sp(2N)))")
    print("The total number of variables is the dimension of this manifold.")
    print("dim(M_b) = dim(SO(4N)/(SO(2N)xSO(2N))) + dim(Sp(4N)/(Sp(2N)xSp(2N)))\n")
    
    # --- Helper functions for dimensions ---
    def dim_so(n):
        return n * (n - 1) // 2

    def dim_sp(m):
        # The group is Sp(m) = USp(m). Dimension is k(2k+1) for USp(2k).
        # So for Sp(m), k = m/2.
        if m % 2 != 0:
            raise ValueError("Argument to Sp must be even.")
        k = m // 2
        return k * (2 * k + 1)
        
    # --- SO Sector Calculation ---
    print("Step 2: Calculate the dimension of the SO sector.")
    so_n_large = 4 * N_replicas
    so_n_small = 2 * N_replicas
    
    dim_so_large = dim_so(so_n_large)
    dim_so_small = dim_so(so_n_small)
    
    print(f"  The SO sector is SO({so_n_large}) / (SO({so_n_small}) x SO({so_n_small})).")
    print(f"  Dimension of SO({so_n_large}) = {so_n_large}({so_n_large}-1)/2 = {dim_so_large}")
    print(f"  Dimension of SO({so_n_small}) = {so_n_small}({so_n_small}-1)/2 = {dim_so_small}")
    
    dim_so_sector = dim_so_large - 2 * dim_so_small
    print(f"  Dimension of SO sector = dim(SO({so_n_large})) - 2 * dim(SO({so_n_small})) = {dim_so_large} - 2 * {dim_so_small} = {dim_so_sector}\n")

    # --- Sp Sector Calculation ---
    print("Step 3: Calculate the dimension of the Sp sector.")
    sp_n_large = 4 * N_replicas
    sp_n_small = 2 * N_replicas

    dim_sp_large = dim_sp(sp_n_large)
    dim_sp_small = dim_sp(sp_n_small)

    print(f"  The Sp sector is Sp({sp_n_large}) / (Sp({sp_n_small}) x Sp({sp_n_small})).")
    print(f"  Dimension of Sp({sp_n_large}) = ({sp_n_large}//2) * (2*({sp_n_large}//2)+1) = {dim_sp_large}")
    print(f"  Dimension of Sp({sp_n_small}) = ({sp_n_small}//2) * (2*({sp_n_small}//2)+1) = {dim_sp_small}")
    
    dim_sp_sector = dim_sp_large - 2 * dim_sp_small
    print(f"  Dimension of Sp sector = dim(Sp({sp_n_large})) - 2 * dim(Sp({sp_n_small})) = {dim_sp_large} - 2 * {dim_sp_small} = {dim_sp_sector}\n")
    
    # --- Total Calculation ---
    print("Step 4: Calculate the total number of variables.")
    total_variables = dim_so_sector + dim_sp_sector
    print("Total non-Grassman variables = (Dimension of SO sector) + (Dimension of Sp sector)")
    print(f"Total variables = {dim_so_sector} + {dim_sp_sector} = {total_variables}")
    
solve_class_d_variables()