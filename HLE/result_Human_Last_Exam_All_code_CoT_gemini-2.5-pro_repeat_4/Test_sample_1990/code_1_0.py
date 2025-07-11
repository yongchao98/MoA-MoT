def efg(q1, q2):
    """
    Simulates the Entanglement Flip Gate (EFG) on two qubits.
    - If states are the same, both flip.
    - If states are different, they remain unchanged.
    """
    if q1 == q2:
        return 1 - q1, 1 - q2
    else:
        return q1, q2

# Initial state of the three-qubit system S = {a, b, c}
# |a> = |0>, |b> = |0>, |c> = |0>
a, b, c = 0, 0, 0
print(f"Initial state: |a, b, c> = |{a}, {b}, {c}>")

# (1) EFG applied to qubits a and b
print("\nStep 1: Applying EFG to a and b")
a, b = efg(a, b)
print(f"State after (1): |a, b, c> = |{a}, {b}, {c}>")

# (2) EFG applied to qubits b and c
print("\nStep 2: Applying EFG to b and c")
b, c = efg(b, c)
print(f"State after (2): |a, b, c> = |{a}, {b}, {c}>")

# (3) EFG applied to qubits a and c
print("\nStep 3: Applying EFG to a and c")
a, c = efg(a, c)
print(f"State after (3): |a, b, c> = |{a}, {b}, {c}>")

# The final state of the system
print("\n-------------------------------------------")
print(f"The final state of the system is |abc> = |{a}{b}{c}>")
print("-------------------------------------------")
<<<|110>>>