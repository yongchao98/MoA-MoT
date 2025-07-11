def efg(q1, q2):
    """
    Applies the Entanglement Flip Gate (EFG) to two qubits.
    - If states are the same, they both flip.
    - If states are different, they remain unchanged.
    States are represented as integers (0 or 1).
    """
    if q1 == q2:
        # Flip to the opposite state (0 becomes 1, 1 becomes 0)
        return 1 - q1, 1 - q2
    else:
        # States remain unchanged
        return q1, q2

# 1. Initial State
a, b, c = 0, 0, 0
print(f"Initial state |abc>: |{a}{b}{c}>")

# 2. EFG applied to qubits a and b
a, b = efg(a, b)
print(f"After EFG on (a,b): |{a}{b}{c}>")

# 3. EFG applied to qubits b and c
b, c = efg(b, c)
print(f"After EFG on (b,c): |{a}{b}{c}>")

# 4. EFG applied to qubits a and c
a, c = efg(a, c)
print(f"After EFG on (a,c): |{a}{b}{c}>")

# Final result
print(f"The final state of the three-qubit system is |{a}{b}{c}>.")
