def efg_gate(q1, q2):
    """Applies the Entanglement Flip Gate to two qubits."""
    if q1 == q2:
        # If states are the same, flip both
        return 1 - q1, 1 - q2
    else:
        # If states are different, they remain unchanged
        return q1, q2

def main():
    """
    Simulates the sequence of EFG operations on a three-qubit system.
    """
    # Initial state: |000>
    qubits = {'a': 0, 'b': 0, 'c': 0}
    
    print(f"Initial state |abc>: |{qubits['a']}{qubits['b']}{qubits['c']}>")
    
    # (1) EFG applied to qubits a and b
    qubits['a'], qubits['b'] = efg_gate(qubits['a'], qubits['b'])
    print(f"State after EFG(a, b): |{qubits['a']}{qubits['b']}{qubits['c']}>")
    
    # (2) EFG applied to qubits b and c
    qubits['b'], qubits['c'] = efg_gate(qubits['b'], qubits['c'])
    print(f"State after EFG(b, c): |{qubits['a']}{qubits['b']}{qubits['c']}>")
    
    # (3) EFG applied to qubits a and c
    qubits['a'], qubits['c'] = efg_gate(qubits['a'], qubits['c'])
    print(f"State after EFG(a, c): |{qubits['a']}{qubits['b']}{qubits['c']}>")
    
    # Final state
    final_state_str = f"|{qubits['a']}{qubits['b']}{qubits['c']}>"
    print(f"\nThe final state of the three-qubit system is: {final_state_str}")

if __name__ == "__main__":
    main()