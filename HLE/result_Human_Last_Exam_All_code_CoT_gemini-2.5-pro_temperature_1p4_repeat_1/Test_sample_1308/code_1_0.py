import numpy as np

# --- Define the Huckel parameters ---
# h_O for Oxygen
h_O = 1.0
# k_CO for Carbon-Oxygen bond
k_CO = 0.8

# --- Construct the Huckel matrix ---
# The molecular structure is O(1)-C(2)-C(3)-O(4)
# The matrix M is set up such that we solve det(M - x*I) = 0
# The eigenvalues of the matrix A below will be our x values.
# A_ii = -h_i (for heteroatoms) or 0 (for carbon)
# A_ij = k_ij (for heteroatom bonds) or 1 (for C-C bonds)
# A_ij = 0 (for non-bonded atoms)
A = np.array([
    [-h_O, k_CO,  0.0,   0.0],
    [k_CO,  0.0,   1.0,   0.0],
    [ 0.0,  1.0,   0.0,  k_CO],
    [ 0.0,  0.0,  k_CO,  -h_O]
])

# --- Calculate eigenvalues ---
# The eigenvalues of matrix A are the 'x' values in E = alpha - x*beta
eigenvalues = np.linalg.eigvals(A)

# --- Sort eigenvalues for energy levels ---
# A larger 'x' value corresponds to a lower, more stable energy.
# We sort them in descending order to list energies from lowest to highest.
x_values = np.sort(eigenvalues)[::-1]

# --- Print the results ---
print("The 4 energies for glyoxal are calculated from the eigenvalues.")
print("The energy E is related to the eigenvalues x by the equation: E = α - xβ")
print("\nThe calculated energies are:")

# Loop through the sorted eigenvalues and print the energy equations
for i, x in enumerate(x_values):
    # Format the output string based on the sign of the eigenvalue x
    if x < 0:
        # If x is negative, E = α - (-|x|)β = α + |x|β
        print(f"E_{i+1} = α + {-x:.3f}β")
    else:
        # If x is positive, E = α - xβ
        print(f"E_{i+1} = α - {x:.3f}β")
