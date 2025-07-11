import numpy as np

# In this script, we model optical components using 2x2 Jones matrices,
# which act on 2x1 Jones vectors representing the light's polarization state.

# --- Case 1: Ideal Birefringent Element (Lossless) ---
# An ideal birefringent medium, like a Quarter-Wave Plate (QWP), is reversible.
# Its Jones matrix is unitary and has an inverse.
print("--- Case 1: Ideal Lossless Birefringent Element (Quarter-Wave Plate) ---")

# Jones matrix for a QWP with its fast axis at 45 degrees.
# It converts linear polarization to circular polarization.
M_qwp = 1/np.sqrt(2) * np.array([[1, -1j],
                                  [-1j, 1]])

print("Jones Matrix for Quarter-Wave Plate (M_qwp):\n", M_qwp)

# A matrix is invertible if and only if its determinant is non-zero.
det_qwp = np.linalg.det(M_qwp)
print(f"\nDeterminant of M_qwp: {det_qwp:.1f}")
print("Since the determinant is non-zero, the matrix is invertible.")

try:
    M_qwp_inv = np.linalg.inv(M_qwp)
    print("\nThe inverse matrix exists, so the operation is reversible in theory.")
    print("Inverse of M_qwp:\n", np.round(M_qwp_inv, 5))
except np.linalg.LinAlgError:
    print("\nCould not compute the inverse.")

print("\n" + "="*70 + "\n")

# --- Case 2: Anisotropic Element with Loss (Polarizer) ---
# A polarizer is an anisotropic element that introduces total loss for one polarization.
# This is an irreversible process.
print("--- Case 2: Element with Irreversible Loss (Linear Polarizer) ---")

# Jones matrix for a perfect horizontal linear polarizer.
# It transmits horizontal light and completely blocks vertical light.
M_pol = np.array([[1.0, 0.0],
                  [0.0, 0.0]])

print("Jones Matrix for Horizontal Polarizer (M_pol):\n", M_pol)

# Check its determinant.
det_pol = np.linalg.det(M_pol)
print(f"\nDeterminant of M_pol: {det_pol}")
print("Since the determinant is zero, the matrix is singular (not invertible).")

print("\nAttempting to calculate the inverse of the polarizer matrix...")
try:
    M_pol_inv = np.linalg.inv(M_pol)
    # This line will not be reached.
    print("Inverse of M_pol:\n", M_pol_inv)
except np.linalg.LinAlgError as e:
    print(f"Result: Failed to compute the inverse. The error is: '{e}'")
    print("\nThe theory fails because the process is not reversible.")
    print("Information about the vertically polarized component of light is permanently lost.")
