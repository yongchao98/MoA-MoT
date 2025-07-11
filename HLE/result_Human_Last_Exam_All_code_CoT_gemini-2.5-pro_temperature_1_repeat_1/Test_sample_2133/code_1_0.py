import numpy as np

def permanent(mat):
    """Calculates the permanent of a 2x2 matrix."""
    if mat.shape != (2, 2):
        raise ValueError("This function only calculates the permanent for 2x2 matrices.")
    return mat[0, 0] * mat[1, 1] + mat[0, 1] * mat[1, 0]

def solve():
    """
    Solves the problem for n=2, which maximizes the specified ratio.
    1. Defines the specific matrix M_2.
    2. Calculates its immanants (determinant and permanent).
    3. Finds and prints the largest immanant.
    """
    # For n=2, the matrix M_n that satisfies the conditions is M_2.
    # It is 2-nilpotent, has all non-zero integer entries, and its
    # Popov form (RREF) maximizes the specified ratio.
    n = 2
    M_n = np.array([[1, 1], 
                    [-1, -1]])

    print(f"The specific matrix for n={n} is M_{n}:")
    print(M_n)
    print("-" * 20)

    # For n=2, the immanants are the determinant and the permanent.
    
    # Calculate the determinant
    det_val = np.linalg.det(M_n)
    print("Calculating the determinant (first immanant):")
    print(f"{M_n[0, 0]} * {M_n[1, 1]} - {M_n[0, 1]} * {M_n[1, 0]} = {det_val}")

    # Calculate the permanent
    perm_val = permanent(M_n)
    print("\nCalculating the permanent (second immanant):")
    print(f"{M_n[0, 0]} * {M_n[1, 1]} + {M_n[0, 1]} * {M_n[1, 0]} = {perm_val}")

    # The largest immanant is the one with the largest absolute value.
    largest_immanant = max(abs(det_val), abs(perm_val))
    print("-" * 20)
    print(f"The largest immanant is the maximum of |{det_val}| and |{perm_val}|.")
    print(f"Result: {largest_immanant}")

solve()