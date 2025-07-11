import sys

def solve_cellular_automaton():
    """
    Finds the intermediate step in a 2-step elementary cellular automaton evolution.
    """
    initial_str = "01101001"
    final_str = "10000111"

    initial_state = [int(c) for c in initial_str]
    final_state = [int(c) for c in final_str]
    length = len(initial_state)

    def apply_rule(state, rule):
        """Applies an ECA rule to a state using periodic boundaries."""
        new_state = []
        state_len = len(state)
        for i in range(state_len):
            # Get neighborhood with periodic boundary conditions
            left = state[(i - 1 + state_len) % state_len]
            center = state[i]
            right = state[(i + 1) % state_len]
            
            # The neighborhood (e.g., 1,0,1) forms a binary number '101' which is 5.
            # This value is used as an index to find the output from the rule's bits.
            # In the standard convention, '111' corresponds to the most significant bit
            # of the rule, and '000' to the least significant.
            neighborhood_value = (left << 2) | (center << 1) | right
            
            # Extract the bit from the rule number that corresponds to the neighborhood
            new_cell_value = (rule >> neighborhood_value) & 1
            new_state.append(new_cell_value)
        return new_state

    # To hold the final answer
    solution_intermediate_state_str = None

    # Iterate through all 256 possible ECA rules
    for rule_num in range(256):
        # First step: initial -> intermediate
        intermediate_candidate = apply_rule(initial_state, rule_num)
        
        # Second step: intermediate -> final
        final_candidate = apply_rule(intermediate_candidate, rule_num)
        
        # Check if the generated final state matches the given final state
        if final_candidate == final_state:
            # Found the unique solution
            solution_intermediate_state_str = "".join(map(str, intermediate_candidate))
            # Since there is only one valid solution, we can stop the search.
            break

    # Print the full sequence as requested by the prompt
    if solution_intermediate_state_str:
        print(initial_str)
        print(solution_intermediate_state_str)
        print(final_str)
    else:
        # This part should not be reached based on the problem description
        print("No solution was found.", file=sys.stderr)
        
solve_cellular_automaton()