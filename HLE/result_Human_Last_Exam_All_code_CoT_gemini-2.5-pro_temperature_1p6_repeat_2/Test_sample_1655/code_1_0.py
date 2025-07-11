import numpy as np

# This script demonstrates why a fixed input polarization and a fixed optical system
# cannot generate an arbitrary output vector polarization state.

# 1. Define the fixed, linearly polarized input beam's polarization.
# We represent horizontal polarization with a Jones vector.
p_in = np.array([[1], [0]])

# 2. Define the transmission matrix T for a single point in the random medium.
# This matrix is a fixed property of the medium. Let's assume some values.
# T = [[T_xx, T_xy], [T_yx, T_yy]]
T = np.array([[0.8 + 0.1j, 0.3 - 0.2j],
              [0.2 + 0.3j, 0.7 - 0.1j]])

# 3. Define the tailored scalar part of the input beam at this point.
# This is our only control. We can change its amplitude and phase.
# Let's just use 1.0 for simplicity.
S_in = 1.0

# 4. Calculate the output electric field vector (unnormalized polarization).
# The equation is: E_out = S_in * T * p_in
E_out = S_in * np.dot(T, p_in)

# 5. Define an "arbitrary" target polarization state we might want to create.
# For example, right-hand circular polarization.
p_target = np.array([[1], [1j]]) / np.sqrt(2) # Normalized target state

# --- Output the results ---
print("--- System Parameters ---")
print(f"Input Polarization Vector (p_in):\n{p_in}")
print("\nTransmission Matrix (T):\n"
      f"[[{T[0,0]:.2f}, {T[0,1]:.2f}],\n"
      f" [{T[1,0]:.2f}, {T[1,1]:.2f}]]")
print(f"\nControlled Scalar Value (S_in): {S_in}")

print("\n--- Calculation ---")
print("Equation: E_out = S_in * T * p_in")
print(f"Resulting Output Field Vector (E_out):\n{E_out}")

# Normalize the output vector to represent the polarization state
p_out = E_out / np.linalg.norm(E_out)
print(f"\nResulting Output Polarization (normalized E_out):\n{p_out}")

print("\n--- Conclusion ---")
print(f"Desired Arbitrary Target Polarization:\n{p_target}")
print("\nNotice that the resulting output polarization is completely determined by the first")
print("column of the matrix T. Changing the scalar S_in would only scale the output")
print("vector, not change its direction (i.e., the polarization state).")
print("Since the resulting polarization is fixed by the system and not equal to our")
print("arbitrary target, we cannot generate any desired vector beam.")
