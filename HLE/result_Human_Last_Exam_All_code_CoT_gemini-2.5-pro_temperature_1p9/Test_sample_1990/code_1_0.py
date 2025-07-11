def efg_gate(q1, q2):
    """
    Applies the Entanglement Flip Gate (EFG) to two qubits.
    - If states are the same, both flip.
    - If states are different, they remain unchanged.
    Returns the new states of the two qubits.
    """
    if q1 == q2:
        # Flip both qubits. 1-q achieves the flip (0->1, 1->0)
        return 1 - q1, 1 - q2
    else:
        # States remain unchanged
        return q1, q2

def main():
    """
    Simulates the sequence of EFG operations on a three-qubit system.
    """
    # Initial state: a=|0>, b=|0>, c=|0>
    a, b, c = 0, 0, 0
    print(f"Initial state of (a, b, c): ({a}, {b}, {c})")

    # (1) EFG applied to qubits a and b
    a, b = efg_gate(a, b)
    print(f"After EFG on (a, b):         ({a}, {b}, {c})")

    # (2) EFG applied to qubits b and c
    b, c = efg_gate(b, c)
    print(f"After EFG on (b, c):         ({a}, {b}, {c})")

    # (3) EFG applied to qubits a and c
    a, c = efg_gate(a, c)
    print(f"After EFG on (a, c):         ({a}, {b}, {c})")

    print("\n--- Final Result ---")
    # As requested, outputting each number in the final state individually
    print(f"Final state of qubit a: {a}")
    print(f"Final state of qubit b: {b}")
    print(f"Final state of qubit c: {c}")

    # Final combined state in ket notation
    final_state_str = f"|{a}{b}{c}>"
    print(f"The final state of the three-qubit system is: {final_state_str}")

if __name__ == "__main__":
    main()
<<<|110>>>