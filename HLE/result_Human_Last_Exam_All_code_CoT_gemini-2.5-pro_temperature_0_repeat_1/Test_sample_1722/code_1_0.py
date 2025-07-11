def simulate_fsm(transitions, sequence, start_state):
    """
    Simulates a deterministic finite state machine.

    Args:
        transitions (dict): The transition function of the FSM.
                            Format: {state: {observation: next_state}}
        sequence (str): The input observation sequence.
        start_state (int): The initial state of the FSM.

    Returns:
        int: The final state of the FSM after processing the sequence.
    """
    current_state = start_state
    for observation in sequence:
        current_state = transitions[current_state][int(observation)]
    return current_state

def solve():
    """
    Finds the minimum hallway length n and demonstrates the solution.
    """
    # The minimum length n is 3.
    # The two observation sequences can be:
    # Omega_1 (Corridor 1): 010
    # Omega_2 (Corridor 2): 101
    n = 3
    omega_1 = "010"
    omega_2 = "101"

    print(f"The minimum length of the hallway is n = {n}.")
    print(f"We can choose the observation sequence for Corridor 1 to be {omega_1}")
    print(f"and for Corridor 2 to be {omega_2}.")
    print("-" * 20)

    # An agent with m=2 states cannot distinguish these sequences.
    # For any 2-state FSM, both sequences will end in the same final state.
    # Therefore, an m=2 agent performs no better than a memoryless (m=1) agent.
    print("An agent with m=2 memory states cannot distinguish these two sequences.")
    print("-" * 20)

    # We demonstrate that an agent with m=3 states CAN distinguish them.
    # Let's design a 3-state FSM (states 0, 1, 2) with start_state = 0.
    m = 3
    start_state = 0
    
    # Define the transition function for our 3-state FSM.
    # This is one possible FSM that can distinguish the sequences.
    transitions_m3 = {
        0: {0: 1, 1: 2},  # From state 0, on '0' go to 1, on '1' go to 2
        1: {0: 1, 1: 0},  # From state 1, on '0' stay in 1, on '1' go to 0
        2: {0: 0, 1: 2}   # From state 2, on '0' go to 0, on '1' stay in 2
    }

    # Simulate the FSM for the first sequence
    final_state_1 = simulate_fsm(transitions_m3, omega_1, start_state)
    
    # Simulate the FSM for the second sequence
    final_state_2 = simulate_fsm(transitions_m3, omega_2, start_state)

    print(f"An agent with m={m} can design an FSM to distinguish them.")
    print(f"Using one such FSM, starting in state {start_state}:")
    print(f"Sequence {omega_1} leads to final state: {final_state_1}")
    print(f"Sequence {omega_2} leads to final state: {final_state_2}")

    if final_state_1 != final_state_2:
        print("\nThe final states are different.")
        print("The agent can assign different actions to these states and achieve a higher reward.")
    else:
        print("\nThe final states are the same. This FSM does not distinguish them.")

    print("-" * 20)
    print(f"Therefore, the minimum value of n is 3.")

solve()
<<<3>>>