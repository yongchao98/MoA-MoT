def entanglement_flip_gate(q1, q2):
    """
    Applies the Entanglement Flip Gate (EFG) operation.
    - If q1 and q2 are the same, they both flip.
    - If q1 and q2 are different, they remain unchanged.
    """
    if q1 == q2:
        # Flip both qubits to the opposite state
        return 1 - q1, 1 - q2
    else:
        # States remain unchanged
        return q1, q2

def solve_qubit_system():
    """
    Solves the state of a three-qubit system after a sequence of EFG applications.
    """
    # Initial state: a=|0>, b=|0>, c=|0>
    a, b, c = 0, 0, 0
    print(f"Initial state (a, b, c): ({a}, {b}, {c})")
    
    # (1) EFG applied to qubits a and b
    print("\nStep 1: Applying EFG to a and b")
    print(f"  - Before: a={a}, b={b}")
    a, b = entanglement_flip_gate(a, b)
    print(f"  - After:  a={a}, b={b}")
    print(f"System state after Step 1: ({a}, {b}, {c})")

    # (2) EFG applied to qubits b and c
    print("\nStep 2: Applying EFG to b and c")
    print(f"  - Before: b={b}, c={c}")
    b, c = entanglement_flip_gate(b, c)
    print(f"  - After:  b={b}, c={c}")
    print(f"System state after Step 2: ({a}, {b}, {c})")
    
    # (3) EFG applied to qubits a and c
    print("\nStep 3: Applying EFG to a and c")
    print(f"  - Before: a={a}, c={c}")
    a, c = entanglement_flip_gate(a, c)
    print(f"  - After:  a={a}, c={c}")
    print(f"System state after Step 3: ({a}, {b}, {c})")

    # Print the final result
    print("\n------------------------------------")
    print(f"The final state of the three-qubit system |a, b, c> is |{a}, {b}, {c}>")
    print("------------------------------------")

# Run the simulation
solve_qubit_system()
