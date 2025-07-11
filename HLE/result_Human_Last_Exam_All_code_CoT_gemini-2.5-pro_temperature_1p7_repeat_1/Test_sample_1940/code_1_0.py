import numpy as np

def print_cccz_construction():
    """
    Prints the construction of a CCCZ gate from 8 CCZ gates.

    This construction is based on an ancilla-free method for building
    a C3-NOT (Toffoli-cubed) gate from 8 C2-NOT (Toffoli) gates.
    The number of CCZ gates required is the same as the number of Toffoli gates.

    Let the qubits be labeled q1, q2, q3 (controls) and q4 (target).
    The gate sequence is presented as an equation.
    """
    
    # Construction from Goudarzi et al. (2020) for a C3-NOT(q1,q2,q3 -> q4)
    # adapted for CCZ gates.
    gate_sequence = [
        "CCZ(q2, q3, q4)",
        "CCZ(q1, q3, q4)",
        "CCZ(q1, q2, q3)",
        "CCZ(q3, q4, q2)",  # Here, q4 acts as a control and q2 as a target
        "CCZ(q1, q2, q3)",
        "CCZ(q3, q4, q2)",
        "CCZ(q1, q3, q4)",
        "CCZ(q2, q3, q4)",
    ]
    
    print("The construction for a CCCZ gate on controls {q1, q2, q3} and target q4 can be expressed as a sequence of 8 CCZ gates:")
    print("\nCCCZ(q1,q2,q3,q4) = ")
    # Print the equation representing the product of gates.
    # Quantum circuits are read from right to left.
    for i, gate in enumerate(reversed(gate_sequence)):
        if i < len(gate_sequence) - 1:
            print(f"  {gate} *")
        else:
            print(f"  {gate}")
            
    num_gates = len(gate_sequence)
    print(f"\nThis construction uses a total of {num_gates} CCZ gates.")
    print("Each number in the final equation represents a qubit index.")
    print("For example, in 'CCZ(q1, q2, q3)', the numbers are 1, 2, and 3.")
    # The construction can be represented by the gate numbers
    # if we map q1->1, q2->2, q3->3, q4->4
    gate_indices = [
        (2, 3, 4),
        (1, 3, 4),
        (1, 2, 3),
        (3, 4, 2),
        (1, 2, 3),
        (3, 4, 2),
        (1, 3, 4),
        (2, 3, 4),
    ]
    
    print("\nRepresenting the qubits by numbers (1, 2, 3, 4):")
    final_equation = " * ".join([f"CCZ{gate}" for gate in reversed(gate_indices)])
    print(f"CCCZ(1,2,3,4) = {final_equation}")

print_cccz_construction()

<<<8>>>