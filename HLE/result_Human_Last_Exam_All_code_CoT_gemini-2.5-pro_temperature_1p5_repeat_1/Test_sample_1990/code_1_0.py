def entanglement_flip_gate(q1, q2):
    """
    Simulates the Entanglement Flip Gate (EFG).
    If qubits are in the same state, they both flip.
    Otherwise, they remain unchanged.
    """
    if q1 == q2:
        return 1 - q1, 1 - q2
    else:
        return q1, q2

def solve_qubit_state():
    """
    Calculates the final state of the three-qubit system.
    """
    # Initial state of the three-qubit system S = {a, b, c}
    a, b, c = 0, 0, 0
    print(f"Initial state:           |{a}{b}{c}>")

    # (1) EFG applied to qubits a and b
    a, b = entanglement_flip_gate(a, b)
    print(f"After EFG on (a, b):     |{a}{b}{c}>")

    # (2) EFG applied to qubits b and c
    b, c = entanglement_flip_gate(b, c)
    print(f"After EFG on (b, c):     |{a}{b}{c}>")

    # (3) EFG applied to qubits a and c
    a, c = entanglement_flip_gate(a, c)
    print(f"Final state (a, b, c):   |{a}{b}{c}>")

if __name__ == "__main__":
    solve_qubit_state()