import math

# --- User-defined parameters ---
# d is the number of sensor nodes (qubits).
# It must be an integer greater than or equal to 1.
d = 1

# F is the fidelity of the noisy GHZ state with respect to the target state.
# It must be a float between 0.0 and 1.0.
F = 0.75

# --- Calculation ---
# The Quantum Fisher Information (QFI) for the parameter theta is given by
# the formula: QFI = 4 * d * (2*F - 1)^2
qfi = 4 * d * (2 * F - 1)**2

# The problem asks for the difference between 1 and the QFI.
result = 1 - qfi

# --- Output ---
# The following print statements show the breakdown of the calculation
# as requested, using the parameters defined above.
print(f"Given parameters: d = {d}, F = {F}")
print(f"The Quantum Fisher Information (QFI) is calculated as: 4 * {d} * (2 * {F} - 1)^2 = {qfi:.4f}")
print(f"The final result (1 - QFI) is: 1 - {qfi:.4f} = {result:.4f}")
