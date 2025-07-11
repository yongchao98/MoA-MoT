def entanglement_flip_gate(q1, q2):
    """
    Applies the Entanglement Flip Gate (EFG) to two qubits.
    - If states are the same, they both flip.
    - If states are different, they remain unchanged.
    """
    if q1 == q2:
        # Flip both qubits (0 -> 1, 1 -> 0)
        return 1 - q1, 1 - q2
    else:
        # States remain unchanged
        return q1, q2

def solve_qubit_state():
    """
    Initializes a three-qubit system and applies a sequence of EFG operations.
    """
    # Initial state: |a, b, c> = |0, 0, 0>
    a, b, c = 0, 0, 0
    print(f"Initial state |a, b, c>: |{a}, {b}, {c}>")

    # (1) EFG applied to qubits a and b
    a, b = entanglement_flip_gate(a, b)
    print(f"After EFG on (a, b): |{a}, {b}, {c}>")

    # (2) EFG applied to qubits b and c
    b, c = entanglement_flip_gate(b, c)
    print(f"After EFG on (b, c): |{a}, {b}, {c}>")

    # (3) EFG applied to qubits a and c
    a, c = entanglement_flip_gate(a, c)
    print(f"Final state |a, b, c>: |{a}, {b}, {c}>")

if __name__ == "__main__":
    solve_qubit_state()