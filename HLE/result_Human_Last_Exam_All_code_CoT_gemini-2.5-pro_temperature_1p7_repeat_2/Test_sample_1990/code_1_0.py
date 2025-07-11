def entanglement_flip_gate(q1, q2):
    """
    Simulates the Entanglement Flip Gate (EFG) on two qubits.

    Args:
        q1 (int): State of the first qubit (0 or 1).
        q2 (int): State of the second qubit (0 or 1).

    Returns:
        tuple: A tuple containing the new states of the two qubits.
    """
    if q1 == q2:
        # w1: if states are the same, flip both
        return (1 - q1, 1 - q2)
    else:
        # w2: if states are different, they remain unchanged
        return (q1, q2)

# Initial state of the three-qubit system S
# a, b, c are all in state |0>
a, b, c = 0, 0, 0
print(f"Initial state: |{a}{b}{c}>")

# (1) EFG applied to qubits a and b
a, b = entanglement_flip_gate(a, b)
print(f"After EFG on a and b: |{a}{b}{c}>")

# (2) EFG applied to qubits b and c
b, c = entanglement_flip_gate(b, c)
print(f"After EFG on b and c: |{a}{b}{c}>")

# (3) EFG applied to qubits a and c
a, c = entanglement_flip_gate(a, c)
print(f"Final state: |{a}{b}{c}>")