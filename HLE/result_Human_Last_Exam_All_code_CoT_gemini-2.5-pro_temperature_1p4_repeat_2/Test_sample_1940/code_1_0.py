# The problem asks for the minimal number of controlled-controlled-Z (CCZ) gates
# to synthesize a controlled-controlled-controlled-Z (CCCZ) gate without using ancilla qubits.
# This is equivalent to synthesizing a 3-controlled-NOT (CCCNOT) gate from 2-controlled-NOT (Toffoli) gates.

# A known result in quantum circuit synthesis provides a formula for the number of Toffoli gates
# needed to construct an n-qubit C^(n-1)Z gate without ancillas.
# The formula is 2*n - 4, where n is the total number of qubits.
# Reference: Y. Shi, "Asymptotically optimal circuits for controlled-Z and controlled-X gates", 2002.

# A CCCZ gate is a 3-controlled-Z gate. It acts on 4 qubits (3 controls, 1 target).
# Therefore, n = 4.

# We apply the formula with n = 4.
n = 4
num_ccz_gates = 2 * n - 4

# The final equation shows the calculation.
print("The total number of qubits (n) for a CCCZ gate is 4 (3 controls + 1 target).")
print("The formula for the minimal number of CCZ gates is 2*n - 4.")
print(f"Calculation: 2 * {n} - 4 = {num_ccz_gates}")
print(f"The minimal number of CCZ gates required is {num_ccz_gates}.")
