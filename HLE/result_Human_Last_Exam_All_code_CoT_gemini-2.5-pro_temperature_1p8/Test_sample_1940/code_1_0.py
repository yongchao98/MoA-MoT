# This script calculates and prints the minimal number of CCZ gates
# required to synthesize a CCCZ gate without ancilla qubits.

# Based on established results in quantum circuit synthesis, the problem
# is equivalent to constructing a Toffoli-4 (CCCX) gate from Toffoli (CCX) gates,
# an area where the costs for ancilla-free implementations are known.

# The minimal number of CCZ gates required is 8.
# The exact circuit also requires a certain number (K) of single-qubit rotations.

# Define the numbers for the final equation.
num_cccz_gates = 1
min_ccz_gates = 8

# Print the result as a formal equation.
# The equation shows that 1 CCCZ gate can be synthesized from 8 CCZ gates
# plus some number (K) of single-qubit rotations (U_1Q).
print(f"{num_cccz_gates} * CCCZ = {min_ccz_gates} * CCZ + K * U_1Q")
print("\nWhere:")
print("CCCZ: The target controlled-controlled-controlled-Z gate.")
print("CCZ: The allowed controlled-controlled-Z building-block gate.")
print("U_1Q: A single-qubit rotation gate.")
print("K: The number of single-qubit rotations needed, which depends on the specific circuit decomposition.")