def solve_automaton():
    """
    Finds the intermediate step in a 2-step elementary cellular automaton evolution
    by testing all possible rules.
    """
    top_row_str = "01101001"
    bottom_row_str = "10000111"

    top_row = [int(c) for c in top_row_str]
    bottom_row = [int(c) for c in bottom_row_str]
    n = len(top_row)
    
    solution_found = False

    # 1. Iterate through all 256 possible rules.
    for rule_num in range(256):
        # The rule number's binary representation defines the automaton's behavior.
        # Per standard convention (e.g., Wolfram's), the binary string b7...b0
        # corresponds to outputs for neighborhoods 111, 110, ..., 000.
        rule_str = format(rule_num, '08b')

        # 2a. Evolve the top row to get a potential middle row.
        potential_mid_row = []
        for i in range(n):
            # Get neighborhood with periodic boundaries.
            left = top_row[(i - 1 + n) % n]
            center = top_row[i]
            right = top_row[(i + 1) % n]
            
            # Convert neighborhood to an integer value (0-7).
            neighborhood_value = left * 4 + center * 2 + right * 1
            
            # Find the new state from the rule string.
            # The bit for neighborhood k is at index (7-k).
            new_state = int(rule_str[7 - neighborhood_value])
            potential_mid_row.append(new_state)

        # 2b. Evolve the potential middle row to get a calculated bottom row.
        calculated_bottom_row = []
        for i in range(n):
            left = potential_mid_row[(i - 1 + n) % n]
            center = potential_mid_row[i]
            right = potential_mid_row[(i + 1) % n]

            neighborhood_value = left * 4 + center * 2 + right * 1
            new_state = int(rule_str[7 - neighborhood_value])
            calculated_bottom_row.append(new_state)
            
        # 3. Check if the result matches the given bottom row.
        if calculated_bottom_row == bottom_row:
            mid_row_str = "".join(map(str, potential_mid_row))
            print(top_row_str)
            print(mid_row_str)
            print(bottom_row_str)
            solution_found = True
            # Since the problem guarantees a unique solution, we can stop.
            return mid_row_str

# Execute the search and print the results.
solution_string = solve_automaton()

# The final answer in the required format.
print(f"<<<{solution_string}>>>")
