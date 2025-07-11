def efg_gate(q1, q2):
    """
    Applies the Entanglement Flip Gate (EFG) to two qubits.
    If the states are the same, they both flip.
    If the states are different, they remain unchanged.
    """
    if q1 == q2:
        return 1 - q1, 1 - q2
    else:
        return q1, q2

def solve():
    """
    Solves the qubit state problem step-by-step.
    """
    # Initial state of the three-qubit system {a, b, c}
    a, b, c = 0, 0, 0
    print(f"Initial state |abc>: |{a}{b}{c}>")

    # (1) EFG applied to qubits a and b
    a, b = efg_gate(a, b)
    print(f"After EFG(a, b), state is: |{a}{b}{c}>")

    # (2) EFG applied to qubits b and c
    b, c = efg_gate(b, c)
    print(f"After EFG(b, c), state is: |{a}{b}{c}>")

    # (3) EFG applied to qubits a and c
    a, c = efg_gate(a, c)
    print(f"After EFG(a, c), state is: |{a}{b}{c}>")

    # Print the final result
    print(f"\nThe final state of the three-qubit system |abc> is |{a}{b}{c}>.")

if __name__ == "__main__":
    solve()
    print("<<<|110>>>>")