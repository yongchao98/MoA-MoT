import numpy as np

# Set a seed for reproducibility of the random matrix
np.random.seed(0)
np.set_printoptions(precision=4, suppress=True)

# 1. Define the initial input beam as horizontally polarized
# We use complex numbers as is standard for Jones Calculus
input_polarization = np.array([[1], [0]], dtype=complex)

# 2. Define the Jones Matrix for the birefringent medium (B)
# Let's model it as a simple 30-degree polarization rotator
theta = np.pi / 6  # 30 degrees
B = np.array([[np.cos(theta), -np.sin(theta)],
              [np.sin(theta),  np.cos(theta)]], dtype=complex)

# 3. Define the Jones Matrix for the random medium (T)
# This is a random 2x2 complex matrix.
T = np.random.rand(2, 2) + 1j * np.random.rand(2, 2)

# Ensure the matrix T is invertible by checking its determinant
if np.linalg.det(T) == 0:
    print("Random matrix T is singular. Please run again.")
    exit()

# Calculate the inverse of the random medium's matrix
T_inv = np.linalg.inv(T)

# --- Theory Test ---

# 4. Propagate the initial input through the system to get Output 1
# The system is B followed by T. The operation is T @ (B @ input_polarization)
output_1 = T @ B @ input_polarization

# 5. Generate Output 2, the proposed new input for the verification step.
# This is done by applying the inverse of the random medium to Output 1.
output_2 = T_inv @ output_1

# 6. Propagate Output 2 through the full system to get the Final Output.
final_output = T @ B @ output_2

# 7. Print and compare the results
# According to the theory, final_output should be equal to output_1.

print("--- Jones Matrices ---")
print("Birefringent Medium (B):\n", B)
print("\nRandom Medium (T):\n", T)

print("\n\n--- Polarization Vectors ---")
print("Initial Input Polarization:\n", input_polarization)
print("\nOutput 1 (after passing Initial Input through T @ B):\n", output_1)
print("\nOutput 2 (the proposed input for the next step, T_inv @ Output_1):\n", output_2)
print("\nFinal Output (after passing Output 2 through T @ B):\n", final_output)

# 8. Check if the theory holds
# We use np.allclose for robust floating-point comparison.
theory_holds = np.allclose(output_1, final_output)

print("\n\n--- Conclusion ---")
print(f"Is Final Output equal to Output 1? \n{theory_holds}")
if not theory_holds:
    print("The theory does not hold. The final output is different from the first output.")
