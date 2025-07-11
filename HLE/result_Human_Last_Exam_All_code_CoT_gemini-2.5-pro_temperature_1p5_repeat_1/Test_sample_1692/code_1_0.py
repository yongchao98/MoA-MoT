def solve_cellular_automaton():
    """
    Finds the intermediate step in an elementary cellular automaton evolution
    by testing all 256 possible rules.
    """
    initial_str = "01101001"
    final_str = "10000111"

    initial_state = [int(c) for c in initial_str]
    final_state = [int(c) for c in final_str]
    n = len(initial_state)

    valid_solutions = []

    for rule_num in range(256):
        # The rule string maps neighborhood 111, 110, ..., 000 to its bits.
        rule_bits = format(rule_num, '08b')

        # --- First generation: initial_state -> intermediate_state ---
        intermediate_state = []
        for i in range(n):
            # Get neighbors with periodic boundary conditions
            left = initial_state[(i - 1 + n) % n]
            center = initial_state[i]
            right = initial_state[(i + 1) % n]

            # Convert neighborhood to an integer (0-7) to find the rule index
            neighborhood_value = 4 * left + 2 * center + 1 * right
            rule_index = 7 - neighborhood_value
            
            new_cell_state = int(rule_bits[rule_index])
            intermediate_state.append(new_cell_state)

        # --- Second generation: intermediate_state -> calculated_final_state ---
        calculated_final_state = []
        for i in range(n):
            # Get neighbors with periodic boundary conditions
            left = intermediate_state[(i - 1 + n) % n]
            center = intermediate_state[i]
            right = intermediate_state[(i + 1) % n]
            
            neighborhood_value = 4 * left + 2 * center + 1 * right
            rule_index = 7 - neighborhood_value
            
            new_cell_state = int(rule_bits[rule_index])
            calculated_final_state.append(new_cell_state)
            
        # Check if the calculated final state matches the target
        if calculated_final_state == final_state:
            valid_solutions.append(intermediate_state)

    # Check for a unique solution and print the result
    if len(valid_solutions) == 0:
        print("No solution was found.")
    else:
        # Remove duplicate solutions if multiple rules produce the same intermediate step
        unique_solutions = []
        for sol in valid_solutions:
            if sol not in unique_solutions:
                unique_solutions.append(sol)
        
        if len(unique_solutions) == 1:
            result_list = unique_solutions[0]
            result_str = "".join(map(str, result_list))
            # The problem requests printing the full equation
            print(f"{initial_str} -> {result_str} -> {final_str}")
        else:
            print("Found multiple possible intermediate steps, so the solution is not unique.")

solve_cellular_automaton()