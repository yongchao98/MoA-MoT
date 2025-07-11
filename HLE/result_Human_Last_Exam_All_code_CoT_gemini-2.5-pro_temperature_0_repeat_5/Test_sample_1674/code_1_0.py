import numpy as np

# --- Setup ---
# 1. Define the initial input beam's polarization (e.g., horizontal)
# This is a Jones vector [Ex, Ey].
input_polarization = np.array([1, 0])

# 2. Define the optical elements as 2x2 Jones matrices.
# For simplicity, we ignore free-space propagation as it doesn't change the core logic.

# A random complex matrix representing the "random medium" T.
np.random.seed(0) # for reproducibility
T = np.random.rand(2, 2) + 1j * np.random.rand(2, 2)
T_inv = np.linalg.inv(T)

# A matrix for a birefringent element B (e.g., a quarter-wave plate).
# This introduces a pi/2 phase shift between the y and x components.
B = np.array([[1, 0], [0, np.exp(1j * np.pi/2)]])

# --- The Experiment ---

# 3. Forward Pass: Propagate the input through the full system (T then B) to get Output 1.
# The light first passes through the random medium T, then the birefringent plate B.
output_1 = B @ T @ input_polarization

# 4. Your Inverse Method: Calculate the theoretical input ("Output 2")
# by applying the inverse of ONLY the random medium (T_inv) to Output 1.
# This is the flawed inverse process described in the prompt.
output_2_calculated = T_inv @ output_1

# 5. Test the Theory: Pass the calculated "Output 2" back through the full system.
# If the theory holds, this result should be identical to Output 1.
reconstructed_output = B @ T @ output_2_calculated

# --- Verification and Output ---

print("This script demonstrates that the proposed theory fails when a birefringent element is added.")
print("-" * 70)

# The final equation we are testing is: reconstructed_output == output_1
# We will print each number in this final comparison.

print(f"Original Output 1 (Ex): {output_1[0]:.4f}")
print(f"Original Output 1 (Ey): {output_1[1]:.4f}")
print("-" * 70)
print(f"Reconstructed Output (Ex): {reconstructed_output[0]:.4f}")
print(f"Reconstructed Output (Ey): {reconstructed_output[1]:.4f}")
print("-" * 70)

# Check if the vectors are close enough to be considered equal.
are_equal = np.allclose(output_1, reconstructed_output)

print(f"Does the reconstructed output match the original Output 1? {are_equal}")
if not are_equal:
    print("The theory does not hold because the inverse process did not account for the birefringent medium.")
