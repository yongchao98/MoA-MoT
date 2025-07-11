def solve_automaton_mystery():
    """
    Finds the intermediate step in a 2-step elementary cellular automaton evolution
    by testing all 256 possible rules.
    """

    def apply_rule(rule_num, state):
        """
        Applies an elementary cellular automaton rule to a given state,
        assuming periodic boundary conditions.
        """
        # The rule's bits correspond to the output for neighborhoods from 111 down to 000.
        # We create a map where the index is the integer value of the neighborhood.
        # E.g., for neighborhood '101' (value 5), the output is the 5th bit of the rule number.
        rule_map = [(rule_num >> i) & 1 for i in range(8)]
        
        n = len(state)
        next_state = [0] * n

        for i in range(n):
            # Get neighborhood with periodic (wraparound) boundaries.
            left = state[(i - 1) % n]
            center = state[i]
            right = state[(i + 1) % n]
            
            # Convert the 3-cell neighborhood pattern to an integer (0-7).
            neighborhood_index = (left << 2) | (center << 1) | right
            
            # Determine the cell's next state from the rule map.
            next_state[i] = rule_map[neighborhood_index]
            
        return next_state

    # The given states in the sequence.
    initial_state_str = "01101001"
    final_state_str = "10000111"
    
    initial_state = [int(c) for c in initial_state_str]
    final_state = [int(c) for c in final_state_str]

    solutions = []

    # Iterate through all 256 possible elementary rules.
    for rule_num in range(256):
        # Step 1: Calculate the intermediate state from the initial state.
        intermediate_state = apply_rule(rule_num, initial_state)
        
        # Step 2: Calculate the final state from the intermediate state.
        calculated_final_state = apply_rule(rule_num, intermediate_state)
        
        # Check if the calculated final state matches the given final state.
        if calculated_final_state == final_state:
            solutions.append(intermediate_state)

    # The problem implies a unique solution, so we print the first one found.
    if solutions:
        solution_state = solutions[0]
        # Print each number of the found intermediate state as a single string.
        print("".join(map(str, solution_state)))

# Execute the solver.
solve_automaton_mystery()