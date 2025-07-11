import math

def calculate_bosonic_variables():
    """
    Calculates the number of non-Grassman variables for the SUSY sigma-model
    in class D with n=2 replicas.
    """
    # The number of replicas is specified as 2.
    n = 2

    # Dimension formula for the orthogonal group O(k).
    def dim_O(k):
        return k * (k - 1) // 2

    # Dimension formula for the compact symplectic group Sp(2k).
    # The formula is dim(Sp(2k)) = k(2k+1).
    def dim_Sp(k_param):
        return k_param * (2 * k_param + 1)

    print("The total number of non-Grassman variables is the dimension of the bosonic target manifold.")
    print("For class D with n replicas, this is dim[O(2n)/(O(n)xO(n))] + dim[Sp(2n)/(Sp(n)xSp(n))].")
    print(f"We are given n = {n}.\n")

    # --- Part 1: Dimension of the compact sector O(4) / (O(2) x O(2)) ---
    print("1. Calculating the dimension of the compact part: O(4) / (O(2) x O(2))")
    
    # Numerator dimension
    k_num_O = 2 * n
    dim_O_num = dim_O(k_num_O)
    print(f"   - Dimension of O({k_num_O}) = {k_num_O}*({k_num_O}-1)/2 = {dim_O_num}")

    # Denominator dimension
    k_den_O = n
    dim_O_den = dim_O(k_den_O)
    dim_O_den_total = 2 * dim_O_den
    print(f"   - Dimension of O({k_den_O}) x O({k_den_O}) = 2 * dim(O({k_den_O})) = 2 * {dim_O_den} = {dim_O_den_total}")

    dim_part1 = dim_O_num - dim_O_den_total
    print(f"   => Dimension of compact part = {dim_O_num} - {dim_O_den_total} = {dim_part1}\n")

    # --- Part 2: Dimension of the non-compact sector Sp(4) / (Sp(2) x Sp(2)) ---
    print("2. Calculating the dimension of the non-compact part: Sp(4) / (Sp(2) x Sp(2))")

    # Numerator dimension
    k_num_Sp = n  # For Sp(2n), the parameter in the formula is n
    dim_Sp_num = dim_Sp(k_num_Sp)
    print(f"   - Dimension of Sp({2*k_num_Sp}) = {k_num_Sp}*(2*{k_num_Sp}+1) = {dim_Sp_num}")

    # Denominator dimension
    k_den_Sp = n // 2 # For Sp(n) where n=2, we have Sp(2), so the parameter is 1.
    dim_Sp_den = dim_Sp(k_den_Sp)
    dim_Sp_den_total = 2 * dim_Sp_den
    print(f"   - Dimension of Sp({2*k_den_Sp}) x Sp({2*k_den_Sp}) = 2 * dim(Sp({2*k_den_Sp})) = 2 * {dim_Sp_den} = {dim_Sp_den_total}")

    dim_part2 = dim_Sp_num - dim_Sp_den_total
    print(f"   => Dimension of non-compact part = {dim_Sp_num} - {dim_Sp_den_total} = {dim_part2}\n")

    # --- Total Dimension ---
    total_dim = dim_part1 + dim_part2
    print("3. Calculating the total number of variables.")
    print("Total variables = (dimension of compact part) + (dimension of non-compact part)")
    print(f"Final Answer = {dim_part1} + {dim_part2} = {total_dim}")

calculate_bosonic_variables()