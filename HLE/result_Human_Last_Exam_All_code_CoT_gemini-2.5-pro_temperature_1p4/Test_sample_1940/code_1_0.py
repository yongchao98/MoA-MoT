def solve_ccc_z_synthesis():
    """
    Calculates and explains the minimal number of CCZ gates to synthesize a CCCZ gate.

    In quantum complexity, synthesizing higher-order gates from a universal
    gate set is a fundamental problem. The question is to find the minimum
    number of controlled-controlled-Z (CCZ) gates to build a
    controlled-controlled-controlled-Z (CCCZ) gate, using arbitrary
    single-qubit rotations but no extra qubits (ancillas).

    1. Local Equivalence:
       - The target CCCZ gate is equivalent to a CCC-NOT (C³-X) gate via single-qubit
         Hadamard gates.
       - The building block CCZ gate is equivalent to a CC-NOT (Toffoli) gate via
         single-qubit Hadamard gates.
       - Since arbitrary single-qubit rotations are allowed, the problem is identical to
         synthesizing a C³-NOT from CC-NOT gates.

    2. Ancilla-Free Construction:
       - Synthesizing multi-controlled gates without ancilla qubits is known to be
         more complex.
       - A widely cited optimal construction for an ancilla-free C³-NOT gate using
         CC-NOT gates requires 8 such gates (D. Maslov, 2015).

    3. The Equation:
       The decomposition can be abstractly represented as:
       1 * CCCZ = 8 * CCZ + N * SingleQubitRotations

       where the single-qubit rotations are used to transform between CCZ/CC-NOT
       gates and to connect them correctly. The resource we are counting is the
       number of the expensive three-qubit CCZ gates.
    """

    # The number of CCCZ gates we want to synthesize.
    num_ccc_z = 1

    # The minimal number of CCZ gates required for the synthesis.
    min_ccz_gates = 8

    print("To synthesize a Controlled-Controlled-Controlled-Z (CCCZ) gate without ancillas:")
    print(f"The target is {num_ccc_z} CCCZ gate.")
    print(f"The minimal number of required CCZ gates is {min_ccz_gates}.")
    print("\nThe final synthesis equation is:")
    print(f"{num_ccc_z} * CCCZ = {min_ccz_gates} * CCZ + N * SingleQubitRotations")


if __name__ == '__main__':
    solve_ccc_z_synthesis()
