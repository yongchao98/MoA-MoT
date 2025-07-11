def solve_qubit_state():
    """
    Simulates the sequence of Entanglement Flip Gate (EFG) operations
    on a three-qubit system to find its final state.
    """
    # Initialize the three qubits a, b, and c to the |0> state.
    # The state is represented by a dictionary.
    qubits = {'a': 0, 'b': 0, 'c': 0}
    print(f"Initial state (a,b,c): |{qubits['a']}{qubits['b']}{qubits['c']}>")
    print("-" * 50)

    def apply_efg(q_name1, q_name2):
        """
        Applies the EFG to two named qubits in the 'qubits' dictionary.
        - q_name1: The name of the first qubit ('a', 'b', or 'c').
        - q_name2: The name of the second qubit ('a', 'b', or 'c').
        """
        q1_val = qubits[q_name1]
        q2_val = qubits[q_name2]
        
        print(f"(Applying EFG to qubits '{q_name1}' and '{q_name2}' with states |{q1_val}> and |{q2_val}> respectively)")
        
        # Rule w1: if states are the same, flip both
        if q1_val == q2_val:
            qubits[q_name1] = 1 - q1_val
            qubits[q_name2] = 1 - q2_val
            print("  -> States were the same. Both qubits flipped.")
        # Rule w2: if states are different, they remain unchanged
        else:
            print("  -> States were different. No change occurred.")
        
        print(f"Current state (a,b,c): |{qubits['a']}{qubits['b']}{qubits['c']}>")
        print("-" * 50)

    # (1) EFG applied to qubits a and b
    apply_efg('a', 'b')

    # (2) EFG applied to qubits b and c
    apply_efg('b', 'c')
    
    # (3) EFG applied to qubits a and c
    apply_efg('a', 'c')

    # Report the final state
    a, b, c = qubits['a'], qubits['b'], qubits['c']
    print("\n--- Simulation Complete ---")
    print(f"The final state of the three-qubit system is |abc> = |{a}{b}{c}>.")
    # The final equation showing each number
    print(f"The final equation is: |a={a}, b={b}, c={c}>")

# Run the simulation
solve_qubit_state()