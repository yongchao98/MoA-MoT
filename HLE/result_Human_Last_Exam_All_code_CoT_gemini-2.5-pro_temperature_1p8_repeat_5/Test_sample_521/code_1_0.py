import sys

def solve_cellular_automaton():
    """
    Analyzes the given cellular automaton pattern to determine its Wolfram code.
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

    # 1. Parse the input string into a grid of 0s and 1s.
    lines = ca_pattern.strip().split('\n')
    grid = [[1 if char == '█' else 0 for char in line.strip()] for line in lines]
    
    height = len(grid)
    width = len(grid[0])
    
    # 2. Deduce the rules from the grid evolution.
    rules = {}
    for r in range(height - 1):
        # Stop once all 8 rules are found.
        if len(rules) == 8:
            break
        for c in range(width):
            # Determine neighborhood with periodic boundary conditions.
            left   = grid[r][(c - 1 + width) % width]
            center = grid[r][c]
            right  = grid[r][(c + 1) % width]
            
            neighborhood = (left, center, right)
            
            # If the rule for this neighborhood is not yet known, record it.
            if neighborhood not in rules:
                result = grid[r + 1][c]
                rules[neighborhood] = result

    # Verify that all 8 rules were found.
    if len(rules) != 8:
        print("Error: Could not determine all 8 rules from the given pattern.", file=sys.stderr)
        return

    # 3. Construct the binary string for the Wolfram code.
    neighborhood_order = [
        (1, 1, 1), (1, 1, 0), (1, 0, 1), (1, 0, 0),
        (0, 1, 1), (0, 1, 0), (0, 0, 1), (0, 0, 0)
    ]
    
    binary_string = "".join(str(rules[n]) for n in neighborhood_order)
        
    # 4. Calculate the decimal value (Wolfram code).
    wolfram_code = int(binary_string, 2)
    
    # 5. Output the results, including the calculation.
    print(f"The binary string representing the rule is: {binary_string}")
    
    calculation_parts = []
    power = len(binary_string) - 1
    for i, digit in enumerate(binary_string):
        if digit == '1':
            value = 2**(power - i)
            calculation_parts.append(str(value))
    
    equation_str = f"{wolfram_code} = " + " + ".join(calculation_parts)
    print(f"The Wolfram code is calculated by converting the binary string to decimal:")
    print(equation_str)

    print(f"\nThus, the Wolfram code is {wolfram_code}.")

solve_cellular_automaton()
<<<75>>>