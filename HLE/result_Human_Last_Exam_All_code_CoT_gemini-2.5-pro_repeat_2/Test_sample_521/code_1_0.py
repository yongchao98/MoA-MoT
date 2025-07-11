import sys

def find_wolfram_code():
    """
    Analyzes a visual representation of a cellular automaton to determine its Wolfram code.
    """
    ca_pattern = """
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
    # 1. Parse the input string into a grid of 0s and 1s
    lines = ca_pattern.strip().split('\n')
    grid = [[1 if char == '█' else 0 for char in line] for line in lines]

    height = len(grid)
    if height < 2:
        print("Error: At least two generations are needed to determine the rule.", file=sys.stderr)
        return

    width = len(grid[0])
    rules = {}

    # 2. Deduce the rules by iterating through the grid
    for r in range(1, height):
        for c in range(width):
            # Assume periodic boundary conditions
            left_parent = grid[r-1][(c - 1 + width) % width]
            center_parent = grid[r-1][c]
            right_parent = grid[r-1][(c + 1) % width]
            
            neighborhood = (left_parent, center_parent, right_parent)
            child_state = grid[r][c]

            # 3. Store the rule if it's new
            if neighborhood not in rules:
                rules[neighborhood] = child_state
            # Optional: Check for inconsistencies
            elif rules[neighborhood] != child_state:
                print(f"Error: Inconsistent rule found for neighborhood {neighborhood}.", file=sys.stderr)
                return
        
        # Optimization: Stop once all 8 rules are found
        if len(rules) == 8:
            break

    if len(rules) != 8:
        print(f"Error: Could not determine all 8 rules from the pattern. Found {len(rules)}.", file=sys.stderr)
        return

    # 4. Construct the Wolfram code
    # Standard order of neighborhoods for Wolfram code (from 111 down to 000)
    neighborhood_order = [
        (1, 1, 1), (1, 1, 0), (1, 0, 1), (1, 0, 0),
        (0, 1, 1), (0, 1, 0), (0, 0, 1), (0, 0, 0)
    ]

    binary_string = "".join([str(rules[n]) for n in neighborhood_order])
    wolfram_code = int(binary_string, 2)

    # 5. Generate the final output showing the calculation
    calculation_parts = []
    for i, bit in enumerate(binary_string):
        power_of_2 = 2**(7-i)
        calculation_parts.append(f"{bit} * {power_of_2}")
    
    calculation_str = " + ".join(calculation_parts)
    print(f"Binary Rule: {binary_string}")
    print(f"Calculation: {calculation_str} = {wolfram_code}")
    print(f"The Wolfram code is: {wolfram_code}")

find_wolfram_code()
<<<77>>>