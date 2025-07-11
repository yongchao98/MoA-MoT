import numpy as np

def solve_cartan_matrix_sum():
    """
    Calculates the sum of entries in the Cartan matrix for the principal block of k[A5 x C2] in char 2.
    """
    # Step 1: Define the decomposition matrix D for the principal block of A5 at p=2.
    # The rows correspond to ordinary characters in the block (dims 1, 3, 3, 5).
    # The columns correspond to the simple modules in the block (dims 1, 2, 2).
    D = np.array([
        [1, 0, 0],
        [1, 1, 0],
        [1, 0, 1],
        [1, 1, 1]
    ])
    print("Step 1: The decomposition matrix D for the principal block of kA5 is:")
    print(D)
    print("-" * 20)

    # Step 2: Calculate the Cartan matrix for the principal block of kA5, C_A5 = D^T * D.
    C_A5 = D.T @ D
    print("Step 2: The Cartan matrix C_A5 for the principal block of kA5 is D^T * D:")
    print(C_A5)
    print("-" * 20)

    # Step 3: The Cartan matrix for the principal block of k[A5 x C2] is 2 * C_A5.
    C_G = 2 * C_A5
    print("Step 3: The Cartan matrix C_G for the principal block of k[A5 x C2] is 2 * C_A5:")
    print(C_G)
    print("-" * 20)

    # Step 4: Calculate the sum of all entries in C_G.
    total_sum = np.sum(C_G)
    
    # Format the output to show the sum equation
    flat_list = C_G.flatten()
    sum_equation = " + ".join(map(str, flat_list))
    
    print("Step 4: The sum of all entries in the Cartan matrix C_G is:")
    print(f"{sum_equation} = {total_sum}")
    
if __name__ == '__main__':
    solve_cartan_matrix_sum()