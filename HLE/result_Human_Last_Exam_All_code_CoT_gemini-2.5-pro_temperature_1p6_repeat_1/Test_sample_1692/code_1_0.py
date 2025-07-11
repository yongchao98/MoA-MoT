def find_intermediate_automaton_step():
    """
    Finds the intermediate step in a 1D elementary cellular automaton sequence
    by searching through all 256 possible rules.
    """
    top_row_str = "01101001"
    bottom_row_str = "10000111"

    # Convert the string representations to lists of integers
    top_row = [int(bit) for bit in top_row_str]
    bottom_row = [int(bit) for bit in bottom_row_str]
    width = len(top_row)

    def apply_rule(generation, rule_number):
        """
        Applies a given ECA rule to a generation to produce the next one.
        Assumes periodic boundary conditions.
        """
        # Get the 8-bit binary representation of the rule
        rule_binary = format(rule_number, '08b')
        next_generation = []

        for i in range(width):
            # Get the left, center, and right neighbors, with wrapping (periodic)
            left = generation[(i - 1 + width) % width]
            center = generation[i]
            right = generation[(i + 1) % width]

            # Convert the 3-cell neighborhood to a number from 0 to 7
            neighborhood_value = 4 * left + 2 * center + 1 * right
            
            # The rule's binary string maps 111, 110, ..., 000 to indices 0, 1, ..., 7
            # So, the new state is determined by the bit at index (7 - neighborhood_value)
            new_state = int(rule_binary[7 - neighborhood_value])
            next_generation.append(new_state)
            
        return next_generation

    # Iterate through all 256 rules to find a match
    solutions = []
    for rule in range(256):
        # Generate the intermediate step from the top row
        intermediate_row = apply_rule(top_row, rule)
        # Generate the final step from the intermediate row
        final_row_candidate = apply_rule(intermediate_row, rule)

        # If the candidate matches the given bottom row, we found a solution
        if final_row_candidate == bottom_row:
            solutions.append("".join(map(str, intermediate_row)))

    # The problem implies a single solution, so we print the first one found.
    if solutions:
        print(solutions[0])
    else:
        print("No solution found.")

find_intermediate_automaton_step()