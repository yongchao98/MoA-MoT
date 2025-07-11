import numpy as np

# Set a consistent print format for complex numbers
np.set_printoptions(precision=3, suppress=True)

print("--- Jones Calculus Demonstration ---")
print("This script simulates how adding a birefringent element breaks the simple inversion theory.\n")

# 1. Define the initial state and optical components using Jones Matrices

# Input: Horizontally polarized light [Ex, Ey]
E_in = np.array([1, 0])
print(f"Initial Input Polarization (E_in):\n{E_in}\n")

# Random Medium T: Assumed to be isotropic (same effect on both polarizations)
# It has some attenuation (0.9) and a phase shift (pi/4 radians)
t_scalar = 0.9 * np.exp(1j * np.pi / 4)
T_matrix = np.array([[t_scalar, 0], [0, t_scalar]])
print(f"Random Medium Matrix (T):\n{T_matrix}\n")

# Birefringent Medium B: A quarter-wave plate at 45 degrees
# This converts horizontal polarization to right-hand circular polarization.
B_matrix = 0.5 * np.array([[1, -1j], [1j, 1]])
print(f"Birefringent Medium Matrix (B) (Quarter-Wave Plate at 45 deg):\n{B_matrix}\n")


# 2. Simulate the forward system to get the target output (Output 1)
# The light passes through the birefringent plate first, then the random medium.
# E_out_1 = T * B * E_in
E_intermediate = B_matrix @ E_in
E_out_1 = T_matrix @ E_intermediate
print(f"--- Calculating Target Output (E_out_1) ---")
print(f"Equation: E_out_1 = T * B * E_in")
print(f"E_out_1:\n{E_out_1}\n")


# 3. Apply the user's inverse theory to find the "necessary input"
# The theory suggests that inverting the random medium's effect on the output
# will yield the necessary input. Let's call this E_supposed_in.
# E_supposed_in = T^-1 * E_out_1

T_inv_matrix = np.linalg.inv(T_matrix)
E_supposed_in = T_inv_matrix @ E_out_1
print(f"--- Applying the Inverse Theory ---")
print(f"Equation: E_supposed_in = T_inv * E_out_1")
print(f"Inverse Random Medium Matrix (T_inv):\n{T_inv_matrix}\n")
print(f"Supposed Necessary Input (E_supposed_in):\n{E_supposed_in}\n")
# Note: E_supposed_in is equal to E_intermediate, which is B * E_in. This makes sense.


# 4. Test the theory: Propagate E_supposed_in through the full system.
# We want to see if this gives us back our target, E_out_1.
# E_final = T * B * E_supposed_in

print(f"--- Testing the Theory ---")
print(f"Equation: E_final = T * B * E_supposed_in")
E_final = T_matrix @ (B_matrix @ E_supposed_in)
print(f"Final calculated output (E_final):\n{E_final}\n")


# 5. Compare the final result with the desired output.
print(f"--- Conclusion ---")
print(f"Target Output (E_out_1):\n{E_out_1}")
print(f"Actual Output (E_final):\n{E_final}")

# Check for equality
are_equal = np.allclose(E_final, E_out_1)
print(f"\nIs the final output equal to the target output? {are_equal}")
if not are_equal:
    print("The theory does not hold. The final output is different from the target.")
    print("This is because the operators for the birefringent medium (B) and the random medium (T) do not commute.")