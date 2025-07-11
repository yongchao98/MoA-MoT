import math

def dim_O(k):
    """Calculates the dimension of the Orthogonal group O(k)."""
    return k * (k - 1) / 2

def dim_Sp(k):
    """Calculates the dimension of the Symplectic group Sp(k). k must be even."""
    if k % 2 != 0:
        raise ValueError("The dimension for Sp(k) must be even.")
    m = k // 2
    return m * (2 * m + 1)

def dim_U(k):
    """Calculates the dimension of the Unitary group U(k)."""
    return k * k

def solve_problem():
    """
    Calculates the number of non-Grassman variables for the supersymmetric sigma-model
    with two replicas for disordered systems of symmetry class D.
    """
    # Number of replicas
    n_r = 2

    # The supergroup G is OSp(2*n_r | 2*n_r) = OSp(4|4)
    # Its bosonic part Gb is O(2*n_r) x Sp(2*n_r) = O(4) x Sp(4)
    dim_O4 = dim_O(4)
    dim_Sp4 = dim_Sp(4)
    dim_Gb = dim_O4 + dim_Sp4

    # The stabilizer subgroup K has a bosonic part Kb = O(n_r) x O(n_r) x U(n_r)
    # For n_r=2, Kb = O(2) x O(2) x U(2)
    dim_O2 = dim_O(2)
    dim_U2 = dim_U(2)
    dim_Kb = dim_O2 + dim_O2 + dim_U2

    # The number of non-Grassman variables is dim(Gb) - dim(Kb)
    num_variables = dim_Gb - dim_Kb
    
    print("The number of non-Grassman variables is the dimension of the bosonic part of the symmetric space G/K.")
    print(f"For n_r = {n_r} replicas in class D, G = OSp({2*n_r}|{2*n_r}) and the bosonic part of the stabilizer K is O({n_r}) x O({n_r}) x U({n_r}).")
    print("\nStep 1: Calculate the dimension of the bosonic part of the supergroup G.")
    print(f"dim(G_b) = dim(O({2*n_r})) + dim(Sp({2*n_r}))")
    print(f"dim(O({2*n_r})) = {int(dim_O4)}")
    print(f"dim(Sp({2*n_r})) = {int(dim_Sp4)}")
    print(f"dim(G_b) = {int(dim_O4)} + {int(dim_Sp4)} = {int(dim_Gb)}")

    print("\nStep 2: Calculate the dimension of the bosonic part of the stabilizer K.")
    print(f"dim(K_b) = dim(O({n_r})) + dim(O({n_r})) + dim(U({n_r}))")
    print(f"dim(O({n_r})) = {int(dim_O2)}")
    print(f"dim(U({n_r})) = {int(dim_U2)}")
    print(f"dim(K_b) = {int(dim_O2)} + {int(dim_O2)} + {int(dim_U2)} = {int(dim_Kb)}")
    
    print("\nStep 3: Calculate the final number of variables.")
    print(f"Number of variables = dim(G_b) - dim(K_b)")
    print(f"Number of variables = {int(dim_Gb)} - {int(dim_Kb)} = {int(num_variables)}")

solve_problem()