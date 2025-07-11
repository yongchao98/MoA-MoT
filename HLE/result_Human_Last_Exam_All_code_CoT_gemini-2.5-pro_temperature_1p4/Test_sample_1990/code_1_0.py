def apply_efg(q1, q2):
    """
    Applies the Entanglement Flip Gate (EFG) to two qubit states.
    - If states are the same (0,0 or 1,1), both flip.
    - If states are different (0,1 or 1,0), they remain unchanged.
    Returns the new states of the two qubits.
    """
    if q1 == q2:
        # Flip both qubits to the opposite state
        return 1 - q1, 1 - q2
    else:
        # Return unchanged states
        return q1, q2

# --- Initialization ---
# The initial state for all qubits {a, b, c} is |0>
a, b, c = 0, 0, 0
print(f"Initial state of the system |abc>: |{a}{b}{c}>")
print("-" * 30)

# --- Sequence of Operations ---

# (1) EFG applied to qubits a and b
print("Step (1): Applying EFG to a and b...")
print(f"  - Before: a={a}, b={b}")
a, b = apply_efg(a, b)
print(f"  - After:  a={a}, b={b}")
print(f"  System state |abc> is now: |{a}{b}{c}>")
print("-" * 30)

# (2) EFG applied to qubits b and c
print("Step (2): Applying EFG to b and c...")
print(f"  - Before: b={b}, c={c}")
b, c = apply_efg(b, c)
print(f"  - After:  b={b}, c={c}")
print(f"  System state |abc> is now: |{a}{b}{c}>")
print("-" * 30)

# (3) EFG applied to qubits a and c
print("Step (3): Applying EFG to a and c...")
print(f"  - Before: a={a}, c={c}")
a, c = apply_efg(a, c)
print(f"  - After:  a={a}, c={c}")
print(f"  System state |abc> is now: |{a}{b}{c}>")
print("-" * 30)

# --- Final Result ---
print(f"After all operations, the final state of the three-qubit system is |{a}{b}{c}>.")

<<<|110>>>