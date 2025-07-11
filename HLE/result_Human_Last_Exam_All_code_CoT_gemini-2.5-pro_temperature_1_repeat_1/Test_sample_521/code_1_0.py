import collections

def solve_cellular_automaton():
    """
    Analyzes the behavior of an elementary cellular automaton to determine its Wolfram code.
    """
    # The visual representation of the automaton's evolution.
    automaton_string = """
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
    lines = automaton_string.strip().split('\n')
    grid = [[1 if char == '█' else 0 for char in line] for line in lines]
    
    height = len(grid)
    width = len(grid[0])
    
    # 2. Extract rules by observing the evolution
    rules = {}
    for r in range(height - 1):
        for c in range(width):
            # Get the 3-cell neighborhood from the current row (r)
            # using periodic (wrap-around) boundary conditions
            left_cell = grid[r][(c - 1 + width) % width]
            center_cell = grid[r][c]
            right_cell = grid[r][(c + 1) + width) % width]
            
            neighborhood = (left_cell, center_cell, right_cell)
            
            # Get the outcome from the next row (r+1)
            outcome = grid[r + 1][c]
            
            # Store the discovered rule
            if neighborhood in rules and rules[neighborhood] != outcome:
                print(f"Error: Inconsistent rule found for neighborhood {neighborhood}.")
                return
            rules[neighborhood] = outcome

    # 3. Build the full rule set in the correct order for the Wolfram code
    neighborhood_order = [
        (1, 1, 1), (1, 1, 0), (1, 0, 1), (1, 0, 0),
        (0, 1, 1), (0, 1, 0), (0, 0, 1), (0, 0, 0)
    ]
    
    if len(rules) != 8:
        print(f"Error: Could not determine all 8 rules. Only found {len(rules)}.")
        print("Missing rules for:")
        for n in neighborhood_order:
            if n not in rules:
                print(n)
        return

    # 4. Calculate the Wolfram code
    binary_rule_str = "".join([str(rules[n]) for n in neighborhood_order])
    wolfram_code = int(binary_rule_str, 2)
    
    # 5. Output the results, including the equation
    print(f"The 8 rules are derived from the pattern.")
    print(f"Ordered binary rule string: {binary_rule_str}")
    
    equation_parts = []
    for i, bit in enumerate(reversed(binary_rule_str)):
        if bit == '1':
            value = 2**i
            equation_parts.append(str(value))
    
    equation = " + ".join(reversed(equation_parts))
    print(f"Decimal Calculation: {equation} = {wolfram_code}")
    print(f"The Wolfram code for this automaton is: {wolfram_code}")

solve_cellular_automaton()