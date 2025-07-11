import numpy as np

# Set up the simulation components

# 1. Initial Input Beam: Horizontally polarized light
# Represented by a Jones Vector: [Ex, Ey]^T
E_in = np.array([[1.0], [0.0]])

# 2. Random Medium (T): A polarization-independent medium.
# It causes some absorption (0.9) and a phase shift (pi/6).
# Since it's isotropic, its matrix is a scalar * identity matrix.
T_scalar = 0.9 * np.exp(1j * np.pi / 6)
T_matrix = T_scalar * np.identity(2)

# Inverse of the Random Medium (T_inv)
T_inv_scalar = 1 / T_scalar
T_inv_matrix = T_inv_scalar * np.identity(2)

# 3. Birefringent Medium (T_b): A half-wave plate at 45 degrees.
# This matrix is NOT a multiple of the identity matrix. It swaps
# horizontal and vertical components, making it strongly polarization-dependent.
T_b_matrix = np.array([[0, 1], [1, 0]], dtype=complex)


# --- Scenario A: Original theory without birefringent medium ---
# Here, we expect the reversal to work perfectly.

print("--- Scenario A: System WITHOUT Birefringent Medium ---")
print(f"Initial Input (E_in):\n{E_in}\n")

# Propagate through the random medium
Output_1A = T_matrix @ E_in
print(f"Output after random medium (T * E_in):\n{Output_1A}\n")

# Reverse the process using the inverse matrix
Output_2A = T_inv_matrix @ Output_1A
print(f"Final State after reversal (T_inv * T * E_in):\n{np.round(Output_2A, 5)}\n")

# Check if the theory holds
print(f"Does the final state equal the initial input? -> {np.allclose(E_in, Output_2A)}")
print("-" * 50)


# --- Scenario B: Add the birefringent medium ---
# Here, we expect the simple reversal (using only T_inv) to fail.

print("\n--- Scenario B: System WITH Birefringent Medium ---")
print(f"Initial Input (E_in):\n{E_in}\n")

# The light first passes through the random medium, then the birefringent plate.
Total_System_Matrix = T_b_matrix @ T_matrix
Output_1B = Total_System_Matrix @ E_in
print(f"Output after combined system (Tb * T * E_in):\n{Output_1B}\n")

# Attempt to reverse the process using ONLY the inverse of the random medium
Output_2B = T_inv_matrix @ Output_1B
print(f"Final State after FAILED reversal (T_inv * Tb * T * E_in):\n{np.round(Output_2B, 5)}\n")

# Check if the theory holds now
print(f"Does the final state equal the initial input? -> {np.allclose(E_in, Output_2B)}")

# As shown by the math, the result should be T_b * E_in, not E_in. Let's verify:
theory_predicts = T_b_matrix @ E_in
print(f"\nThe theory predicts the final state should be Tb * E_in:\n{theory_predicts}")
print(f"Is our result consistent with this prediction? -> {np.allclose(Output_2B, theory_predicts)}")
