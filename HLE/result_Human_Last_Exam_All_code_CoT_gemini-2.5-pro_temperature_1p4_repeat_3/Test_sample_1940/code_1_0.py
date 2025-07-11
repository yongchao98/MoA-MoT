import sys

def solve_cccz_synthesis():
    """
    Calculates the minimal number of CCZ gates to synthesize a CCCZ gate without ancillas.

    This problem is a known case in quantum circuit synthesis. The number of required
    gates can be determined by a general formula for decomposing multi-controlled gates.

    The synthesis of a C^n-Z gate (a Z gate with n controls) from C^(n-1)-Z gates has the
    same gate count as synthesizing a C^n-NOT from C^(n-1)-NOT gates, because one can be
    converted to the other via single-qubit Hadamard gates on the target.

    For ancilla-free synthesis, a C^n-NOT gate can be constructed from 2n-2
    C^(n-1)-NOT gates. This construction is known to be optimal.

    In this problem, we are synthesizing a CCCZ (C^3-Z) gate, so n=3.
    The building blocks are CCZ (C^2-Z) gates.
    """

    # The number of controls on the target gate (CCCZ) is 3.
    n = 3

    # The formula for the minimal number of C^(n-1) gates is 2n-2.
    num_ccz_gates = 2 * n - 2

    # The problem asks us to output the numbers in the final equation.
    # The equation represents the synthesis:
    # 1 CCCZ = (minimal number of CCZ gates) * CCZ + (single-qubit rotations for corrections)

    print("This problem asks for the minimal number of CCZ gates to exactly synthesize one CCCZ gate without using ancilla qubits.")
    print("Based on established results in quantum circuit synthesis, we can use a general formula.\n")
    print(f"The target CCCZ gate is a C^n-Z gate with n = {n}.")
    print(f"The minimal number of C^(n-1)-Z gates (CCZ gates) required for an ancilla-free synthesis is given by the formula 2*n - 2.")
    print(f"Calculation: 2 * {n} - 2 = {num_ccz_gates}\n")
    print("This leads to the following synthesis equation:")
    
    cccz_count = 1
    print(f"{cccz_count} * CCCZ = {num_ccz_gates} * CCZ + single-qubit rotations")

solve_cccz_synthesis()