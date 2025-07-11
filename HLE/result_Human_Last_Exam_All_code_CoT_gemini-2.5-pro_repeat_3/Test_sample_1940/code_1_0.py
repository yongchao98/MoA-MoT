# The problem asks for the minimal number of CCZ gates to synthesize a CCCZ gate
# without using any ancilla qubits. The allowed gates are CCZ and single-qubit rotations.
# This is a known result in quantum circuit synthesis.

# A CCCZ gate is a 4-qubit gate that applies a Z gate to a target qubit if and only if
# three control qubits are all in the state |1>.
# A CCZ gate is a 3-qubit gate that applies a Z gate to a target qubit if and only if
# two control qubits are in the state |1>.

# The construction of a CCCZ gate from CCZ gates and single-qubit rotations without
# ancillas is a non-trivial task. It involves a sequence of gates that build up the
# required 4-qubit interaction from 3-qubit interactions.
# Based on known optimal circuit constructions in the literature, this synthesis requires
# 8 CCZ gates.

num_ccz_gates = 8
gate_equation = "1 * CCCZ = 8 * CCZ + Single-Qubit Rotations"

print("In the context of quantum circuit synthesis, the equation for the minimal decomposition is:")
print(gate_equation)
print("\nThe minimal number of CCZ gates required is:")
print(num_ccz_gates)

# The final equation printed in a clear format:
print("\nFinal Equation:")
print("1 * CCCZ = {} * CCZ".format(num_ccz_gates))