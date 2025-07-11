import numpy as np

# Set print options for clarity
np.set_printoptions(precision=3, suppress=True)

# 1. Define the input beam (purely horizontal polarization)
E_in = np.array([[1 + 0j], [0 + 0j]])
print(f"Input Polarization Vector (E_in):\n{E_in}\n")

# 2. Define a simple random medium 'T' that mixes polarizations.
# For realism, we make it symmetric (reciprocal)
# T_ij represents coupling from polarization j to polarization i
T = np.array([
    [0.8+0.1j, 0.3-0.2j],
    [0.3-0.2j, 0.7-0.1j]
])
print(f"Random Medium Matrix (T):\n{T}\n")

# 3. Define a birefringent medium 'J_b' with different phase shifts
# for horizontal (phi_x) and vertical (phi_y) components.
phi_x = np.pi / 4  # 45 degree phase shift
phi_y = np.pi / 2  # 90 degree phase shift
J_b = np.array([
    [np.exp(1j*phi_x), 0],
    [0, np.exp(1j*phi_y)]
])
print(f"Birefringent Medium Matrix (J_b):\n{J_b}\n")

# 4. Case 1: Light passes through the random medium T, then birefringence J_b.
# Final Output = J_b * T * E_in
E_out_1 = J_b @ T @ E_in

print("--- Case 1: T followed by J_b ---")
# Final Equation: E_out_1 = J_b * T * E_in
print(f"Equation: E_out_1 = \n{J_b}\n * \n{T}\n * \n{E_in}")
print(f"\nFinal Output Vector (E_out_1):\n{E_out_1}\n")


# 5. Case 2: Light passes through birefringence J_b, then the random medium T.
# Final Output = T * J_b * E_in
E_out_2 = T @ J_b @ E_in

print("--- Case 2: J_b followed by T ---")
# Final Equation: E_out_2 = T * J_b * E_in
print(f"Equation: E_out_2 = \n{T}\n * \n{J_b}\n * \n{E_in}")
print(f"\nFinal Output Vector (E_out_2):\n{E_out_2}\n")


# 6. Conclusion
# We use np.allclose to check if the two output vectors are equal.
are_equal = np.allclose(E_out_1, E_out_2)
print("--- Conclusion ---")
print(f"Are the outputs from Case 1 and Case 2 equal? {are_equal}")
print("Since the final outputs are different, the order of operations matters.")
print("This means T and J_b do not commute, breaking any theory that assumes the system's inverse is unchanged.")
