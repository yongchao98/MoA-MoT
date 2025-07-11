import math

# --- Step 1: Define the parameters of the problem ---
# d: the total number of sensor nodes (qubits).
#    Let's use an example value, e.g., d = 5.
d = 5

# F: the fidelity of the noisy GHZ state with respect to the pure GHZ state
#    (|0...0> + |1...1>)/sqrt(2). F must be between 0 and 1.
#    Let's use an example value, e.g., F = 0.9.
F = 0.9

print(f"Calculating for d = {d} and F = {F}\n")

# --- Step 2: Calculate the Quantum Fisher Information (QFI) ---
# The derived formula for the QFI is: I_Q = 4 * d * (2*F - 1)^2
# Let's break down the calculation as requested.

# The equation for the final result is: 1 - QFI
# where QFI = term1 * term2 * term3
# term1 = 4
# term2 = d
# term3 = (2*F - 1)^2

term1 = 4
term2 = d
term3_base = 2 * F - 1
term3 = term3_base**2

# Calculate the QFI
qfi = term1 * term2 * term3

# --- Step 3: Calculate the final result ---
# The problem asks for the difference between 1 and the QFI.
final_result = 1 - qfi

# --- Step 4: Print the components of the final equation ---
print("The final result is calculated from the equation: 1 - QFI")
print("where QFI = 4 * d * (2*F - 1)^2\n")

print("--- Components of the equation ---")
print(f"The constant term '1': {1}")
print(f"The constant '4' in the QFI formula: {term1}")
print(f"The number of sensors 'd': {term2}")
print(f"The fidelity-dependent term '(2*F - 1)^2': ({2}*{F} - 1)^2 = {term3_base:.2f}^2 = {term3:.4f}")
print("-------------------------------------\n")

print(f"Calculated Quantum Fisher Information (QFI): {qfi:.4f}")
print(f"Final result (1 - QFI): {final_result:.4f}")
