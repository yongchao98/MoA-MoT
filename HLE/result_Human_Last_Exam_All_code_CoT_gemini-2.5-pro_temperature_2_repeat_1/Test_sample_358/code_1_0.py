import numpy as np

def solve_cartan_sum():
    """
    Calculates the sum of entries in the Cartan matrix of the principal block of k[A5 x C2] for p=2.
    """
    print("The group is G = A_5 x C_2, and the characteristic of the field is 2.")
    print("The sum of entries in the Cartan matrix of the principal block of G is the product of the sums for A_5 and C_2.")
    print("-" * 30)

    # Step 1: Cartan matrix for the principal block of A5 at p=2.
    # The principal block B0(A5) has ordinary characters chi_1, chi_2, chi_3, chi_5
    # and simple modules S1 (dim 1) and S2 (dim 2).
    # The decomposition matrix D, expressing ordinary characters in terms of simple (Brauer) characters, is:
    D_A5 = np.array([
        [1, 0],  # chi_1 -> 1*S1 + 0*S2
        [1, 1],  # chi_2 -> 1*S1 + 1*S2
        [1, 1],  # chi_3 -> 1*S1 + 1*S2
        [1, 2]   # chi_5 -> 1*S1 + 2*S2
    ])
    
    # The Cartan matrix C is given by C = D^T * D.
    C0_A5 = D_A5.T @ D_A5
    sum_A5 = np.sum(C0_A5)
    
    print("1. For A_5:")
    print("The Cartan matrix of the principal block is:")
    print(str(C0_A5))
    c_a5_11 = C0_A5[0, 0]
    c_a5_12 = C0_A5[0, 1]
    c_a5_21 = C0_A5[1, 0]
    c_a5_22 = C0_A5[1, 1]
    print(f"The sum of entries for A_5 is {c_a5_11} + {c_a5_12} + {c_a5_21} + {c_a5_22} = {int(sum_A5)}.")
    print("-" * 30)
    
    # Step 2: Cartan matrix for the principal block of C2 at p=2.
    # For G = C2, the group algebra kC2 is the principal block. There is one simple module S0.
    # The projective cover P0 is kC2 itself, which has composition factors S0, S0.
    # Thus, the Cartan matrix is the 1x1 matrix C0(C2) = [[2]].
    C0_C2 = np.array([[2]])
    sum_C2 = np.sum(C0_C2)

    print("2. For C_2:")
    c_c2_11 = C0_C2[0,0]
    print(f"The Cartan matrix of the principal block is [[{c_c2_11}]].")
    print(f"The sum of entries for C_2 is {int(sum_C2)}.")
    print("-" * 30)
    
    # Step 3: Combine the results.
    total_sum = sum_A5 * sum_C2
    
    print("3. Final Calculation:")
    print(f"The total sum of entries is the product of the individual sums.")
    print(f"Total Sum = ({c_a5_11} + {c_a5_12} + {c_a5_21} + {c_a5_22}) * ({c_c2_11}) = {int(sum_A5)} * {int(sum_C2)} = {int(total_sum)}")

solve_cartan_sum()