def efg(q1, q2):
    """
    Applies the Entanglement Flip Gate (EFG) operation on two qubits.
    - If states are the same, they both flip.
    - If states are different, they remain unchanged.
    """
    if q1 == q2:
        return (1 - q1, 1 - q2)
    else:
        return (q1, q2)

def print_state(step_name, a, b, c):
    """Prints the current state of the qubits."""
    print(f"{step_name}: a={a}, b={b}, c={c}")

# 1. Initial state
a, b, c = 0, 0, 0
print_state("Initial State          ", a, b, c)

# 2. First operation: EFG applied to qubits a and b
a, b = efg(a, b)
print_state("After EFG(a, b)        ", a, b, c)

# 3. Second operation: EFG applied to qubits b and c
b, c = efg(b, c)
print_state("After EFG(b, c)        ", a, b, c)

# 4. Third operation: EFG applied to qubits a and c
a, c = efg(a, c)
print_state("Final State (after EFG(a, c))", a, b, c)

print("\nThe final state of the three-qubit system (a,b,c) is |" + str(a) + str(b) + str(c) + ">")
