def solve_qubit_state():
    """
    Calculates the final state of a three-qubit system after a series of
    Entanglement Flip Gate (EFG) operations.
    """
    # Initial state of the three-qubit system {a, b, c} is |000>
    a, b, c = 0, 0, 0
    print(f"Initial state: |a,b,c> = |{a},{b},{c}>")

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

    # (1) EFG applied to qubits a and b
    # a=0, b=0 are the same, so they flip to 1, 1
    a, b = efg(a, b)
    print(f"After EFG(a, b): |a,b,c> = |{a},{b},{c}>")

    # (2) EFG applied to qubits b and c
    # b=1, c=0 are different, so they remain unchanged
    b, c = efg(b, c)
    print(f"After EFG(b, c): |a,b,c> = |{a},{b},{c}>")

    # (3) EFG applied to qubits a and c
    # a=1, c=0 are different, so they remain unchanged
    a, c = efg(a, c)
    print(f"After EFG(a, c): |a,b,c> = |{a},{b},{c}>")
    
    print("\nThe final state of the three-qubit system is:")
    # The final equation requested is the representation of the final state |abc>
    print(f"|{a}{b}{c}>")

solve_qubit_state()