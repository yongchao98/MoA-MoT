import math

def calculate_variables():
    """
    Calculates the number of non-Grassman variables for the supersymmetric sigma-model
    with two replicas for disordered systems of symmetry class D.
    """
    # Number of replicas
    n = 2
    print(f"The number of replicas is n = {n}.\n")

    # The target manifold is M = G/H where G = O(2n|2n) and H = O(n|n) x O(n|n).
    # The number of bosonic variables is dim_B(M) = dim_B(G) - dim_B(H).

    # Define a helper function for the dimension of an orthogonal group O(m)
    def dim_O(m):
        return m * (m - 1) // 2

    # --- Step 1: Calculate the dimension of the bosonic part of G = O(2n|2n) ---
    p_G = 2 * n
    q_G = 2 * n
    print(f"The supergroup is G = O({p_G}|{q_G}).")
    
    # dim_B(G) = dim(O(2n)) + dim(O(2n))
    dim_O_pG = dim_O(p_G)
    dim_O_qG = dim_O(q_G)
    dim_b_G = dim_O_pG + dim_O_qG
    
    print(f"The dimension of its bosonic part is dim_B(G) = dim(O({p_G})) + dim(O({q_G})) = {dim_O_pG} + {dim_O_qG} = {dim_b_G}.\n")

    # --- Step 2: Calculate the dimension of the bosonic part of H = O(n|n) x O(n|n) ---
    p_H_factor = n
    q_H_factor = n
    print(f"The subgroup is H = O({p_H_factor}|{q_H_factor}) x O({p_H_factor}|{q_H_factor}).")

    # First, find the dimension for a single O(n|n) factor
    # dim_B(O(n|n)) = dim(O(n)) + dim(O(n))
    dim_O_pH_factor = dim_O(p_H_factor)
    dim_O_qH_factor = dim_O(q_H_factor)
    dim_b_H_factor = dim_O_pH_factor + dim_O_qH_factor
    
    # The total dimension for H is twice the dimension of one factor
    dim_b_H = 2 * dim_b_H_factor
    
    print(f"The dimension of the bosonic part of one O({p_H_factor}|{q_H_factor}) factor is dim(O({p_H_factor})) + dim(O({q_H_factor})) = {dim_O_pH_factor} + {dim_O_qH_factor} = {dim_b_H_factor}.")
    print(f"The total dimension of the bosonic part of H is 2 * {dim_b_H_factor} = {dim_b_H}.\n")

    # --- Step 3: Calculate the final result ---
    num_variables = dim_b_G - dim_b_H
    
    print("The number of non-Grassman (bosonic) variables is the difference:")
    print(f"Result = dim_B(G) - dim_B(H) = {dim_b_G} - {dim_b_H} = {num_variables}")

if __name__ == "__main__":
    calculate_variables()