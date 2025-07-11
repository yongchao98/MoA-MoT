import numpy as np

# Set print options for clarity
np.set_printoptions(precision=3, suppress=True)

print("--- Step 1: Define the optical components and input beam ---\n")

# Input Beam: Horizontally polarized light [Ex, Ey]
E_in = np.array([[1], [0]], dtype=complex)
print(f"Initial Input Beam (E_in):\n{E_in}\n")

# Original Medium T: A simple, non-polarizing phase plate (adds a 45-degree phase shift)
# It acts the same on both horizontal and vertical polarization.
phase_shift = np.exp(1j * np.pi / 4)
T = np.array([[phase_shift, 0], [0, phase_shift]], dtype=complex)
T_inv = np.linalg.inv(T)
print(f"Original Medium Matrix (T):\n{T}\n")
print(f"Inverse of Original Medium (T^-1):\n{T_inv}\n")

# Birefringent Medium B: A quarter-wave plate at 45 degrees
# This medium affects different polarizations differently.
B = 0.5 * np.array([[1-1j, 1+1j], [1+1j, 1-1j]], dtype=complex)
print(f"Birefringent Medium Matrix (B):\n{B}\n")

# The new system includes both media. The new transmission matrix T_new is their product.
T_new = B @ T
print(f"New Combined Medium Matrix (T_new = B @ T):\n{T_new}\n")


print("\n--- Step 2: Propagate through the new system to get Output 1 ---\n")
# Output_1 = T_new * E_in
Output_1 = T_new @ E_in
print(f"The final output of the new system is (Output_1 = T_new @ E_in):\n{Output_1}\n")


print("\n--- Step 3: Apply the original theory's inversion step ---\n")
# The theory suggests using the inverse of the original medium (T^-1) to find the required input.
# Predicted_Input = T^-1 * Output_1
Predicted_Input = T_inv @ Output_1
print("The theory predicts that the necessary input to produce Output_1 is:")
print(f"(Predicted_Input = T^-1 @ Output_1):\n{Predicted_Input}\n")
print("Note: This is not the same as the original E_in.\n")


print("\n--- Step 4: Test the theory by using the Predicted_Input ---\n")
# If the theory holds, propagating Predicted_Input through the new system (T_new)
# should result in the original Output_1.
# Final_Output = T_new * Predicted_Input
Final_Output = T_new @ Predicted_Input

print("We test this by calculating: Final_Output = T_new @ Predicted_Input")
print("Equation components:")
print(f"T_new:\n{T_new}\n")
print(f"Predicted_Input:\n{Predicted_Input}\n")
print(f"Result (Final_Output):\n{Final_Output}\n")


print("\n--- Step 5: Conclusion ---\n")
# Compare the Final_Output with the original Output_1
print(f"Original Output_1 was:\n{Output_1}\n")
print(f"Final_Output from predicted input is:\n{Final_Output}\n")

# The np.allclose function checks if two arrays are element-wise equal within a tolerance.
are_equal = np.allclose(Final_Output, Output_1)

print(f"Does the Final_Output match the original Output_1? -> {are_equal}\n")
if not are_equal:
    print("Conclusion: The theory does not hold. The calculated Final_Output is different from the actual Output_1.")
    print("This is because the inversion step used T^-1, which does not account for the birefringent medium B.")

<<<Yes, the theory can fail. The birefringent medium introduces a polarization-dependent effect that is not captured by the inverse of the original transmission matrix. Therefore, the inversion step becomes incorrect for the new system.>>>