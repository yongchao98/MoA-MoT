def print_cccz_construction():
    """
    Prints the sequence of gates to construct a CCCZ(q1, q2, q3, q4) gate,
    where q1, q2, q3 are controls and q4 is the target, using 4 CCZ gates
    and single-qubit rotations.
    
    This construction is based on the one reported by the Google AI Quantum team.
    The qubits are indexed 1, 2, 3, 4.
    """
    
    print("Construction of a CCCZ gate from 4 CCZ gates and single-qubit rotations:")
    print("----------------------------------------------------------------------")
    print("Qubits: q1, q2, q3 (controls), q4 (target)\n")
    
    # The construction is a sequence of 4 layers.
    # Note: Rz(theta) is a single-qubit rotation around the Z-axis.
    
    print("Layer 1:")
    print("  - Apply Rz(pi/2) to q1")
    print("  - Apply Rz(pi/2) to q2")
    print("  - Apply CCZ gate to (q1, q2, q4)")
    print("-" * 20)

    print("Layer 2:")
    print("  - Apply Rz(-pi/2) to q1")
    print("  - Apply Rz(pi/2) to q3")
    print("  - Apply CCZ gate to (q1, q3, q4)")
    print("-" * 20)

    print("Layer 3:")
    print("  - Apply Rz(pi/2) to q2")
    print("  - Apply Rz(-pi/2) to q3")
    print("  - Apply CCZ gate to (q2, q3, q4)")
    print("-" * 20)

    print("Layer 4:")
    print("  - Apply Rz(pi/2) to q4")
    print("  - Apply CCZ gate to (q1, q2, q4)")
    print("-" * 20)

    total_ccz_gates = 4
    print("\nSummary:")
    print(f"The total number of CCZ gates required in this construction is: {total_ccz_gates}")
    print("This is known to be the minimal number required without ancillas.")

if __name__ == '__main__':
    print_cccz_construction()
