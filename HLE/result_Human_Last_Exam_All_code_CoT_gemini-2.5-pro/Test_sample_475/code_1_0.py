import numpy as np
from scipy.special import lambertw

# Define the constants from the problem
sigma_0 = 7.43e-7  # units: e/nm
R_0 = 30.0         # units: nm
q_i = 2 * np.pi

# --- Calculation Step 1: Compute the constant C ---
# C = W(exp(q_i)) / (1 + W(exp(q_i)))^3
# The Lambert W function lambertw() can return complex numbers,
# but for a positive real argument exp(q_i), the principal branch (k=0) is real.
exp_qi = np.exp(q_i)
W_exp_qi = lambertw(exp_qi).real

# Ensure W_exp_qi is not -1 to avoid division by zero, which is not an issue here.
if 1 + W_exp_qi == 0:
    raise ValueError("Denominator is zero in constant C calculation.")

C = W_exp_qi / ((1 + W_exp_qi)**3)

# --- Calculation Step 2: Compute pi^3 ---
pi_cubed = np.pi**3

# --- Calculation Step 3: Compute the total charge Q ---
# Formula: Q = sigma_0 * R_0 * pi^3 * C
Q = sigma_0 * R_0 * pi_cubed * C

# --- Output the results ---
# As requested, output each number in the final equation
print("The final calculation is based on the formula: Q = sigma_0 * R_0 * pi^3 * C")
print("where C = W(exp(q_i)) / (1 + W(exp(q_i)))^3")
print("\n--- Values used in the equation ---")
print(f"sigma_0 (e/nm) = {sigma_0}")
print(f"R_0 (nm) = {R_0}")
print(f"q_i = {q_i}")
print(f"pi^3 = {pi_cubed}")
print(f"C = {C}")

# Print the final equation with all values substituted
print("\n--- Final Equation ---")
print(f"Q = {sigma_0} * {R_0} * {pi_cubed} * {C}")

# Print the final answer for the total charge
print("\n--- Total Charge ---")
print(f"The total charge on the droplet is Q = {Q:.6e} e")

# Final answer in the required format
# Let's round the answer to 2 decimal places in scientific notation
answer_str = f"{Q:.2e}"
# Converting to a number format, e.g. 1.76e-05 -> 1.76E-5 -> 0.0000176
final_answer_val = float(answer_str)
# Let's check the magnitude. 1.76e-5 seems reasonable. Let's output it as a decimal.
# print(f"<<<{final_answer_val}>>>")
# The question format is <<<answer content>>>. The number itself is fine.
print(f"\n<<<%.5e>>>" % Q)