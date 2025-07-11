def efg(qubits, q1_key, q2_key):
    """
    Applies the Entanglement Flip Gate (EFG) to two qubits in the system.
    
    Args:
    qubits (dict): The dictionary representing the state of the qubit system.
    q1_key (str): The key for the first qubit.
    q2_key (str): The key for the second qubit.
    """
    print(f"Applying EFG to qubits {q1_key} and {q2_key}.")
    print(f"  - State before: a={qubits['a']}, b={qubits['b']}, c={qubits['c']}")
    
    # Check if the states of the two qubits are the same
    if qubits[q1_key] == qubits[q2_key]:
        print(f"  - Qubits {q1_key} ({qubits[q1_key]}) and {q2_key} ({qubits[q2_key]}) are in the same state. Flipping both.")
        # Flip both qubits. new_state = 1 - current_state
        qubits[q1_key] = 1 - qubits[q1_key]
        qubits[q2_key] = 1 - qubits[q2_key]
    else:
        print(f"  - Qubits {q1_key} ({qubits[q1_key]}) and {q2_key} ({qubits[q2_key]}) are in different states. No change.")
        
    print(f"  - State after:  a={qubits['a']}, b={qubits['b']}, c={qubits['c']}\n")

def solve_qubit_puzzle():
    """
    Solves the qubit puzzle by simulating the sequence of EFG operations.
    """
    # Initial state: |abc> = |000>
    qubits = {'a': 0, 'b': 0, 'c': 0}
    print(f"Initial state: a={qubits['a']}, b={qubits['b']}, c={qubits['c']}\n")

    # Sequence of EFG applications
    # (1) EFG applied to qubits a and b
    efg(qubits, 'a', 'b')
    
    # (2) EFG applied to qubits b and c
    efg(qubits, 'b', 'c')
    
    # (3) EFG applied to qubits a and c
    efg(qubits, 'a', 'c')

    # Output the final state
    print("Final state of the three-qubit system:")
    final_state_str = f"|{qubits['a']}{qubits['b']}{qubits['c']}>"
    print(f"a = {qubits['a']}")
    print(f"b = {qubits['b']}")
    print(f"c = {qubits['c']}")
    print(f"Which corresponds to the state: {final_state_str}")

# Run the simulation
solve_qubit_puzzle()