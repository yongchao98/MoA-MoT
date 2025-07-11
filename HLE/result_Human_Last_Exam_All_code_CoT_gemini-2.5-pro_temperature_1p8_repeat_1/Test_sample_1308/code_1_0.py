import numpy as np

def solve_glyoxal_huckel():
    """
    Calculates the 4 pi-system energies of glyoxal using modified Huckel theory.

    The problem provides the following parameters for the heteroatom (Oxygen):
    h_O = 1
    k_CO = 0.8

    For carbon atoms, the standard parameters are used:
    h_C = 0
    k_CC = 1.0

    The energy levels E are given by E = α + λβ, where λ are the eigenvalues of
    the simplified Huckel matrix B.
    """
    # 1. Define the Huckel parameters from the problem description.
    h_O = 1.0
    k_CO = 0.8
    h_C = 0.0
    k_CC = 1.0

    # 2. Construct the simplified 4x4 Huckel matrix B for O(1)-C(2)-C(3)-O(4).
    # B[i, j] = h_i if i == j
    # B[i, j] = k_ij if atoms i and j are bonded
    # B[i, j] = 0 otherwise
    B = np.array([
        [h_O,  k_CO, 0.0,  0.0 ],  # Row for O(1)
        [k_CO, h_C,  k_CC, 0.0 ],  # Row for C(2)
        [0.0,  k_CC, h_C,  k_CO],  # Row for C(3)
        [0.0,  0.0,  k_CO, h_O ]   # Row for O(4)
    ])

    # 3. Calculate the eigenvalues of the matrix B.
    # These eigenvalues are the 'λ' values.
    eigenvalues = np.linalg.eigvals(B)

    # 4. Sort the eigenvalues in ascending order.
    eigenvalues_sorted = np.sort(eigenvalues)

    # 5. Print the final results in the specified format.
    print("The 4 energies (E) of glyoxal are given by the formula E = α + λβ.")
    print("The calculated values for λ are:")
    for i, val in enumerate(eigenvalues_sorted):
        # Format the output to show the final equation clearly
        if val < 0:
            print(f"E{i+1} = α - {abs(val):.4f}β")
        else:
            print(f"E{i+1} = α + {val:.4f}β")

# Execute the function to find and print the energies.
solve_glyoxal_huckel()