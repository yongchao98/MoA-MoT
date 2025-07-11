def entanglement_flip_gate(q1, q2):
    """
    Applies the Entanglement Flip Gate (EFG) to two qubits.
    - If states are the same, they both flip.
    - If states are different, they remain unchanged.
    """
    if q1 == q2:
        return 1 - q1, 1 - q2
    else:
        return q1, q2

# Initial state of the three-qubit system {a, b, c}
a, b, c = 0, 0, 0

print(f"Initial state (a, b, c) = ({a}, {b}, {c})")
print("-" * 30)

# (1) EFG applied to qubits a and b
a, b = entanglement_flip_gate(a, b)
print(f"State after EFG(a, b):  (a, b, c) = ({a}, {b}, {c})")

# (2) EFG applied to qubits b and c
b, c = entanglement_flip_gate(b, c)
print(f"State after EFG(b, c):  (a, b, c) = ({a}, {b}, {c})")

# (3) EFG applied to qubits a and c
a, c = entanglement_flip_gate(a, c)
print(f"State after EFG(a, c):  (a, b, c) = ({a}, {b}, {c})")

print("-" * 30)
# The final state is represented by the integer values of the qubits.
# In quantum notation, this corresponds to the state |abc>.
print(f"The final state of the three-qubit system is (a, b, c) = ({a}, {b}, {c})")
print(f"In quantum notation, the final state is |{a}{b}{c}>")
