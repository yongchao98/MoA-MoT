def apply_efg(q1, q2):
    """
    Applies the Entanglement Flip Gate (EFG) to two qubits.
    - If the qubits are in the same state, they both flip.
    - If they are in different states, they remain unchanged.
    Returns the new states of the two qubits.
    """
    if q1 == q2:
        return 1 - q1, 1 - q2
    else:
        return q1, q2

# Initial state of the three-qubit system S = {a, b, c}
a, b, c = 0, 0, 0
print(f"Initial state |abc>: |{a}{b}{c}>")

# (1) EFG applied to qubits a and b
# a=0, b=0 -> a=1, b=1
a, b = apply_efg(a, b)
print(f"After EFG(a, b), state |abc> is: |{a}{b}{c}>")

# (2) EFG applied to qubits b and c
# b=1, c=0 -> b=1, c=0 (unchanged)
b, c = apply_efg(b, c)
print(f"After EFG(b, c), state |abc> is: |{a}{b}{c}>")

# (3) EFG applied to qubits a and c
# a=1, c=0 -> a=1, c=0 (unchanged)
a, c = apply_efg(a, c)
print(f"After EFG(a, c), state |abc> is: |{a}{b}{c}>")

# Final result
print(f"\nThe final state of the three-qubit system is |{a}{b}{c}>.")

<<<|110>>>