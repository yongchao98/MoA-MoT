def entanglement_flip_gate(q1, q2):
    """
    Applies the Entanglement Flip Gate (EFG) to two qubits.
    - If states are the same, both flip.
    - If states are different, they remain unchanged.
    Returns the new states of the two qubits.
    """
    if q1 == q2:
        return 1 - q1, 1 - q2
    else:
        return q1, q2

# Initial state of the three-qubit system S = {a, b, c}
qubits = {'a': 0, 'b': 0, 'c': 0}

print(f"Initial state: a={qubits['a']}, b={qubits['b']}, c={qubits['c']}")

# (1) EFG applied to qubits a and b
print("\n(1) Applying EFG to a and b...")
qubits['a'], qubits['b'] = entanglement_flip_gate(qubits['a'], qubits['b'])
print(f"State after (1): a={qubits['a']}, b={qubits['b']}, c={qubits['c']}")

# (2) EFG applied to qubits b and c
print("\n(2) Applying EFG to b and c...")
qubits['b'], qubits['c'] = entanglement_flip_gate(qubits['b'], qubits['c'])
print(f"State after (2): a={qubits['a']}, b={qubits['b']}, c={qubits['c']}")

# (3) EFG applied to qubits a and c
print("\n(3) Applying EFG to a and c...")
qubits['a'], qubits['c'] = entanglement_flip_gate(qubits['a'], qubits['c'])
print(f"State after (3): a={qubits['a']}, b={qubits['b']}, c={qubits['c']}")

# Print the final state
print("\n--------------------------------")
print("Final state of the three-qubit system {a, b, c}:")
final_a = qubits['a']
final_b = qubits['b']
final_c = qubits['c']
print(f"|{final_a}{final_b}{final_c}>")
print("--------------------------------")

<<<110>>>