def solve_cellular_automaton():
    """
    Finds the unique intermediate step in a 3-generation elementary cellular automaton sequence.
    """
    top_row = '01101001'
    target_bottom_row = '10000111'
    
    def apply_rule(rule_binary, input_row):
        """
        Applies a given ECA rule to an input row with cyclic boundaries.
        """
        n = len(input_row)
        next_row = []
        
        # A mapping from the neighborhood (e.g., '111') to the rule's output bit.
        # The rule_binary string corresponds to outputs for neighborhoods 
        # '111', '110', '101', '100', '011', '010', '001', '000' in that order.
        rule_map = {f'{i:03b}': rule_binary[7 - i] for i in range(8)}

        for i in range(n):
            # Get neighbors with cyclic boundary conditions
            left = input_row[(i - 1 + n) % n]
            center = input_row[i]
            right = input_row[(i + 1) % n]
            
            neighborhood = left + center + right
            next_row.append(rule_map[neighborhood])
            
        return "".join(next_row)

    valid_solutions = []
    
    # Iterate through all 256 possible rules
    for rule_num in range(256):
        # Format the rule number as an 8-bit binary string
        rule_binary = format(rule_num, '08b')
        
        # 1. Generate the candidate middle row from the top row
        middle_row_candidate = apply_rule(rule_binary, top_row)
        
        # 2. Generate a test bottom row from the candidate middle row
        generated_bottom_row = apply_rule(rule_binary, middle_row_candidate)
        
        # 3. Check if the generated bottom row matches the target
        if generated_bottom_row == target_bottom_row:
            valid_solutions.append(middle_row_candidate)
            
    # Find the unique solution among the valid ones found
    unique_solutions = list(set(valid_solutions))
    
    if len(unique_solutions) == 1:
        solution = unique_solutions[0]
        print("The given sequence is:")
        print(f"Top row:    {top_row}")
        print(f"Middle row: ????????")
        print(f"Bottom row: {target_bottom_row}")
        print("\nThe only valid solution for the intermediate step is:")
        # The problem asks to "output each number in the final equation"
        # We will print the solution string, which is composed of numbers.
        print(solution)
    elif len(unique_solutions) == 0:
        print("No valid solution was found.")
    else:
        print("Multiple possible solutions were found:")
        for sol in unique_solutions:
            print(sol)

solve_cellular_automaton()
<<<10101010>>>