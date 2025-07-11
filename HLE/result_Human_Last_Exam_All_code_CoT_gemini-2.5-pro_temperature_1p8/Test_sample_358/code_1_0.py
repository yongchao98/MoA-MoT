import numpy as np

def solve_cartan_sum():
    """
    Calculates the sum of all entries in the Cartan matrix of the principal block
    of the group algebra kG, where G = A_5 x C_2 and char(k) = 2.
    """

    # --- Introduction ---
    print("This program calculates the sum of the Cartan matrix entries for the principal block of k[A_5 x C_2].")
    print("The strategy is to compute the sums for the principal blocks of kC_2 and kA_5 separately and then multiply them.")
    print("-" * 30)

    # --- Part 1: Principal Block of kC_2 ---
    print("Step 1: Analyzing the principal block of kC_2 for p=2.")
    # For char(k)=2, kC_2 has one block (the principal one) with one simple module.
    # The algebra is local, and its Cartan matrix is [[2]].
    cartan_C2 = np.array([[2]])
    sum_C2 = np.sum(cartan_C2)
    print("The Cartan matrix for the principal block of kC_2 is:")
    print(cartan_C2)
    print(f"The sum of its entries is: {sum_C2}")
    print("-" * 30)

    # --- Part 2: Principal Block of kA_5 ---
    print("Step 2: Analyzing the principal block of kA_5 for p=2.")
    # From the representation theory of A_5, we know the decomposition matrix (D) for the principal 2-block.
    # The rows correspond to the ordinary characters in the block (chi_1, chi_2, chi_3, chi_5).
    # The columns correspond to the simple modules in the block (S_1, S_2, S_3).
    D_A5 = np.array([
        [1, 0, 0],  # chi_1 restricts to S_1
        [1, 1, 0],  # chi_2 restricts to S_1 + S_2
        [1, 0, 1],  # chi_3 restricts to S_1 + S_3
        [1, 1, 1]   # chi_5 restricts to S_1 + S_2 + S_3
    ])
    print("The decomposition matrix D for the principal 2-block of kA_5 is:")
    print(D_A5)
    
    # The Cartan matrix C is calculated as C = D^T * D.
    cartan_A5 = D_A5.T @ D_A5
    print("\nThe Cartan matrix C = D^T * D for the principal 2-block of kA_5 is:")
    print(cartan_A5)

    # Summing the entries of the Cartan matrix for A_5.
    sum_A5 = np.sum(cartan_A5)
    print(f"\nThe sum of the entries of the Cartan matrix for kA_5 is: {sum_A5}")
    print("-" * 30)
    
    # --- Part 3: Final Calculation ---
    print("Step 3: Calculating the total sum for k[A_5 x C_2].")
    total_sum = sum_A5 * sum_C2
    print("The total sum is the product of the individual sums.")
    print(f"Total Sum = (Sum for kA_5) * (Sum for kC_2)")
    print(f"            = {sum_A5} * {sum_C2} = {total_sum}")

solve_cartan_sum()