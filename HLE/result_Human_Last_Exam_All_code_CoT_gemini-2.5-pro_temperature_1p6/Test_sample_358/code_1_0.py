import numpy as np

def solve():
    """
    Calculates the sum of all entries in the Cartan matrix of the principal block
    of the group algebra kG, where G = A5 x C2 and char(k) = 2.
    """

    # Step 1: Define the decomposition matrix for the principal block of A5 at p=2.
    # The rows correspond to ordinary characters in the block (dims 1, 3, 3, 5).
    # The columns correspond to simple modules in the block (dims 1, 2, 2).
    D_A5 = np.array([
        [1, 0, 0],  # chi_1 = phi_1
        [1, 1, 0],  # chi_2 = phi_1 + phi_2
        [1, 0, 1],  # chi_3 = phi_1 + phi_3
        [1, 1, 1]   # chi_5 = phi_1 + phi_2 + phi_3
    ])

    # Step 2: Compute the Cartan matrix for the principal block of A5.
    # C_A5 = D_A5^T * D_A5
    C_A5 = D_A5.transpose() @ D_A5

    print("The Cartan matrix for the principal block of A5 at p=2 is:")
    print(C_A5)

    # Step 3: Sum the entries of the Cartan matrix for A5.
    sum_A5 = np.sum(C_A5)
    
    print("\nThe sum of its entries is:")
    C_A5_flat = C_A5.flatten()
    print(" + ".join(map(str, C_A5_flat)), f"= {int(sum_A5)}")


    # Step 4: Define the Cartan matrix for the principal block of C2 at p=2.
    # kC2 is uniserial with composition factors [S, S], where S is the trivial module.
    C_C2 = np.array([[2]])

    print("\nThe Cartan matrix for the principal block of C2 at p=2 is:")
    print(C_C2)

    # Step 5: Sum the entries of the Cartan matrix for C2.
    sum_C2 = np.sum(C_C2)
    print("\nThe sum of its entries is:")
    print(f"{int(sum_C2)} = {int(sum_C2)}")

    # Step 6: The total sum for G = A5 x C2 is the product of the individual sums.
    total_sum = sum_A5 * sum_C2

    print("\nThe total sum of entries for the Cartan matrix of the principal block of k[A5 x C2] is the product of these sums:")
    print(f"{int(sum_A5)} * {int(sum_C2)} = {int(total_sum)}")

solve()