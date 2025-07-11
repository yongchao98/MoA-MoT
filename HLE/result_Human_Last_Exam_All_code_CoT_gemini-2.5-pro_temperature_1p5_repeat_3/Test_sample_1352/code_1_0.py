import itertools
from collections import Counter

def solve_switch_problem():
    """
    Calculates the expected number of rounds for the switches to return to their initial state.
    """
    # Step 1: Define influence sets based on the problem description.
    # Person number (1-8) mapped to the set of people they influence.
    influence_sets = {
        1: {2, 4, 6, 7}, 2: {3, 5, 6, 8}, 3: {4, 6},
        4: {5},           5: {6, 8},       6: {7},
        7: {8},           8: {}
    }

    # Convert influence sets to 0-indexed influence vectors for easier computation.
    influence_vectors = {}
    for person, targets in influence_sets.items():
        vec = [0] * 8
        for target in targets:
            vec[target - 1] = 1
        influence_vectors[person] = tuple(vec)

    # Memoization cache to speed up the transition function.
    transition_cache = {}

    # Step 2: Define the function for a single round's state transition.
    def transition_function(state):
        """
        Calculates the next state after one round, given the current state.
        State is a tuple of 8 integers (0 or 1).
        """
        if state in transition_cache:
            return transition_cache[state]

        new_state_list = list(state)
        # People take turns in sequence from person 8 down to 1.
        for person in range(8, 0, -1):
            person_idx = person - 1
            # If the person's switch is ON at the time of their turn...
            if new_state_list[person_idx] == 1:
                # ...all people in their influence set flip their switches.
                v_i = influence_vectors[person]
                for j in range(8):
                    new_state_list[j] = (new_state_list[j] + v_i[j]) % 2
        
        result = tuple(new_state_list)
        transition_cache[state] = result
        return result

    # Step 3: Find all cycles and their lengths.
    total_states = 2**8
    all_states = list(itertools.product([0, 1], repeat=8))
    
    visited_states = set()
    cycle_lengths = []
    
    for state in all_states:
        if state in visited_states:
            continue

        # Found a new cycle. Trace it to find its length.
        cycle = []
        current_state = state
        # The path must lead back to the start state, forming a cycle.
        while current_state not in visited_states:
            visited_states.add(current_state)
            cycle.append(current_state)
            current_state = transition_function(current_state)
        
        cycle_lengths.append(len(cycle))

    # Step 4: Calculate the expected value E[R].
    # For a cycle of length r, all r states in it have a return time of r.
    # The total sum of return times is the sum of r^2 for each cycle.
    total_return_time_sum = sum(r * r for r in cycle_lengths)
    
    expected_value = total_return_time_sum / total_states
    
    # Step 5: Format and print the final output.
    cycle_counts = Counter(cycle_lengths)
    
    print("The state space is partitioned into the following cycles:")
    for length, count in sorted(cycle_counts.items()):
        print(f"- {count} cycle(s) of length {length}")
        
    print("\nThe expected value E[R] is the average of the return times over all 256 states.")
    print("The total sum of return times is calculated by summing r*r for each cycle of length r.")
    
    equation_parts = []
    sum_parts = []
    for length, count in sorted(cycle_counts.items()):
        term_value = (length ** 2) * count
        equation_parts.append(f"({length}^2 * {count})")
        sum_parts.append(str(term_value))
    
    numerator_calc = " + ".join(equation_parts)
    numerator_sum = " + ".join(sum_parts)
    
    print(f"\nE[R] = ( {numerator_calc} ) / {total_states}")
    if len(sum_parts) > 1:
        print(f"     = ( {numerator_sum} ) / {total_states}")
    print(f"     = {total_return_time_sum} / {total_states}")
    print(f"     = {expected_value}")
    
    print(f"\nRounded to 2 decimal places, the expected value is: {expected_value:.2f}")


solve_switch_problem()
<<<29.75>>>