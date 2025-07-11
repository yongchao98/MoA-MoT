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
qubits = {'a': 0, 'b': 0, 'c': 0}
print(f"Initial state |abc>: |{qubits['a']}{qubits['b']}{qubits['c']}>")
print("-" * 30)

# (1) EFG applied to qubits a and b
print("1. Applying EFG to a and b...")
qubits['a'], qubits['b'] = apply_efg(qubits['a'], qubits['b'])
print(f"   State after step (1): |{qubits['a']}{qubits['b']}{qubits['c']}>")
print("-" * 30)

# (2) EFG applied to qubits b and c
print("2. Applying EFG to b and c...")
qubits['b'], qubits['c'] = apply_efg(qubits['b'], qubits['c'])
print(f"   State after step (2): |{qubits['a']}{qubits['b']}{qubits['c']}>")
print("-" * 30)

# (3) EFG applied to qubits a and c
print("3. Applying EFG to a and c...")
qubits['a'], qubits['c'] = apply_efg(qubits['a'], qubits['c'])
print(f"   State after step (3): |{qubits['a']}{qubits['b']}{qubits['c']}>")
print("-" * 30)

# Print the final state
print("After all operations, the final state of the three-qubit system is:")
# Printing each number in the final equation as requested
a, b, c = qubits['a'], qubits['b'], qubits['c']
print(f"|{a}{b}{c}>")
