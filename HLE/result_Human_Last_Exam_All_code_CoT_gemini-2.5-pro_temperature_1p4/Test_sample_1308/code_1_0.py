import numpy as np

def calculate_huckel_energies_glyoxal():
    """
    Calculates the 4 pi-electron energies of glyoxal using modified Huckel theory.
    
    The secular determinant is solved as an eigenvalue problem. The energies are
    expressed in the form E = α + xβ, where x are the eigenvalues of the
    Huckel matrix B.
    
    For glyoxal (O1-C2-C3-O4), with parameters h_O=1, h_C=0, k_CO=0.8, k_CC=1,
    the matrix B is:
    [[h_O, k_CO, 0,   0  ]
     [k_CO, h_C, k_CC, 0  ]
     [0,   k_CC, h_C, k_CO]
     [0,   0,   k_CO, h_O ]]
    """
    # Heteroatom parameters
    h_O = 1.0
    k_O = 0.8
    h_C = 0.0
    k_CC = 1.0

    # Construct the Huckel matrix B
    B = np.array([
        [h_O, k_O, 0.0, 0.0],
        [k_O, h_C, k_CC, 0.0],
        [0.0, k_CC, h_C, k_O],
        [0.0, 0.0, k_O, h_O]
    ])

    # The values of x in E = α + xβ are the eigenvalues of B
    eigenvalues = np.linalg.eigvals(B)

    # Sort eigenvalues in descending order. Since β is a negative energy unit,
    # the largest x value corresponds to the lowest energy level (most stable).
    sorted_eigenvalues = np.sort(eigenvalues)[::-1]

    print("The four π-electron energies of glyoxal are (from lowest to highest):")
    
    # E1 = α + (x1)β
    x1 = sorted_eigenvalues[0]
    print(f"E1 = α + {x1:.4f}β")

    # E2 = α + (x2)β
    x2 = sorted_eigenvalues[1]
    print(f"E2 = α + {x2:.4f}β")

    # E3 = α + (x3)β, where x3 is negative
    x3 = sorted_eigenvalues[2]
    print(f"E3 = α - {-x3:.4f}β")
    
    # E4 = α + (x4)β, where x4 is negative
    x4 = sorted_eigenvalues[3]
    print(f"E4 = α - {-x4:.4f}β")

if __name__ == "__main__":
    calculate_huckel_energies_glyoxal()