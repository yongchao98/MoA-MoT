import math

def evaluate_integral(zeta_i, zeta_j):
  """
  Evaluates the integral <phi_i | 1/r | phi_j> for 1s Slater orbitals
  based on the analytical formula.

  Args:
    zeta_i: The orbital exponent for the first 1s Slater orbital (phi_i).
    zeta_j: The orbital exponent for the second 1s Slater orbital (phi_j).

  Returns:
    The value of the integral.
  """
  numerator = 4 * (zeta_i * zeta_j)**(1.5)
  denominator = (zeta_i + zeta_j)**2
  result = numerator / denominator
  return result

# --- Main execution ---
# This part of the script will be executed when run.

# Case 1: Diagonal element where i = j.
# This corresponds to the expectation value of the 1/r operator.
# Let's use zeta = 1.0, the optimal value for a hydrogen atom.
zeta_1 = 1.0
print("Evaluating the diagonal element <phi | 1/r | phi>:")
print(f"For a single 1s orbital with zeta = {zeta_1}:")
result_1 = evaluate_integral(zeta_1, zeta_1)
# The output demonstrates how the final number is calculated
print(f"Result from formula: 4 * ({zeta_1} * {zeta_1})**1.5 / ({zeta_1} + {zeta_1})**2 = {result_1}")
print(f"As shown in the derivation, the result simplifies to zeta, which is {zeta_1}.")
print("-" * 30)

# Case 2: Off-diagonal element where i != j.
# This integral appears in calculations using a basis set of multiple STOs.
# Let's use two different zeta values as an example.
zeta_i = 1.0
zeta_j = 1.2
print("Evaluating the off-diagonal element <phi_i | 1/r | phi_j>:")
print(f"For two 1s orbitals with zeta_i = {zeta_i} and zeta_j = {zeta_j}:")
result_2 = evaluate_integral(zeta_i, zeta_j)
# The output demonstrates how the final number is calculated
print(f"Result from formula: 4 * ({zeta_i} * {zeta_j})**1.5 / ({zeta_i} + {zeta_j})**2 = {result_2:.6f}")
