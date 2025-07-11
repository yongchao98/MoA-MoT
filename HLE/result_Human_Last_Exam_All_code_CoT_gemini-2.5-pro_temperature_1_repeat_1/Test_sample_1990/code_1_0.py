def entanglement_flip_gate(q1, q2):
    """
    Applies the Entanglement Flip Gate (EFG) to two qubits.
    - If states are the same, they both flip.
    - If states are different, they remain unchanged.
    """
    if q1 == q2:
        return 1 - q1, 1 - q2
    else:
        return q1, q2

def solve_qubit_state():
    """
    Initializes a three-qubit system and applies a sequence of EFG operations.
    """
    # Initial state of the three-qubit system S = {a, b, c}
    a, b, c = 0, 0, 0
    print(f"Initial state (a, b, c): |{a}{b}{c}>")

    # (1) EFG applied to qubits a and b
    print("\n(1) Applying EFG to a and b...")
    a, b = entanglement_flip_gate(a, b)
    print(f"State after (1): |{a}{b}{c}>")

    # (2) EFG applied to qubits b and c
    print("\n(2) Applying EFG to b and c...")
    b, c = entanglement_flip_gate(b, c)
    print(f"State after (2): |{a}{b}{c}>")

    # (3) EFG applied to qubits a and c
    print("\n(3) Applying EFG to a and c...")
    a, c = entanglement_flip_gate(a, c)
    print(f"State after (3): |{a}{b}{c}>")

    # Final result
    print(f"\nThe final state of the three-qubit system is for a, b, and c is {a}, {b}, and {c} respectively.")
    print(f"This corresponds to the state |{a}{b}{c}>.")

if __name__ == "__main__":
    solve_qubit_state()