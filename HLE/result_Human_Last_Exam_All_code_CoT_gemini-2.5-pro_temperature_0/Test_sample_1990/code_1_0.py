def efg(q1, q2):
    """
    Applies the Entanglement Flip Gate (EFG) to two qubits.
    - If states are the same, they both flip.
    - If states are different, they remain unchanged.
    """
    if q1 == q2:
        return 1 - q1, 1 - q2
    else:
        return q1, q2

def main():
    """
    Simulates the sequence of EFG applications on a three-qubit system.
    """
    # Initial state of the three-qubit system S = {a, b, c}
    a, b, c = 0, 0, 0
    print(f"Initial state |a, b, c>: |{a}, {b}, {c}>")

    # (1) EFG applied to qubits a and b
    a, b = efg(a, b)
    print(f"After EFG(a, b), state is: |{a}, {b}, {c}>")

    # (2) EFG applied to qubits b and c
    b, c = efg(b, c)
    print(f"After EFG(b, c), state is: |{a}, {b}, {c}>")

    # (3) EFG applied to qubits a and c
    a, c = efg(a, c)
    print(f"After EFG(a, c), state is: |{a}, {b}, {c}>")

    print("\nFinal state of the three-qubit system:")
    print(f"|a, b, c> = |{a}, {b}, {c}>")

if __name__ == "__main__":
    main()
<<<|1, 1, 0>>>