def apply_efg(q1, q2):
    """
    Applies the Entanglement Flip Gate (EFG) to two qubits.
    - If states are the same, they both flip.
    - If states are different, they remain unchanged.
    """
    if q1 == q2:
        return 1 - q1, 1 - q2
    else:
        return q1, q2

# Initial state of the three-qubit system S = {a, b, c}
a, b, c = 0, 0, 0
print(f"Initial state: |a,b,c> = |{a},{b},{c}>")

# (1) EFG applied to qubits a and b
print("\nApplying EFG to qubits a and b...")
a, b = apply_efg(a, b)
print(f"State after (1): |a,b,c> = |{a},{b},{c}>")

# (2) EFG applied to qubits b and c
print("\nApplying EFG to qubits b and c...")
b, c = apply_efg(b, c)
print(f"State after (2): |a,b,c> = |{a},{b},{c}>")

# (3) EFG applied to qubits a and c
print("\nApplying EFG to qubits a and c...")
a, c = apply_efg(a, c)
print(f"State after (3): |a,b,c> = |{a},{b},{c}>")

# Final result
print(f"\nAfter all operations, the final state of the three-qubit system is |{a}{b}{c}>.")
print(f"This is composed of the final states for each qubit: a={a}, b={b}, c={c}")
<<<|110>>>