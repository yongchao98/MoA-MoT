import textwrap

def solve_cellular_automaton():
    """
    Determines the Wolfram code for an elementary cellular automaton from its visual evolution.
    """
    automaton_pattern = """
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
    lines = textwrap.dedent(automaton_pattern).strip().split('\n')
    grid = [[1 if char == '█' else 0 for char in line] for line in lines]
    
    height = len(grid)
    if height < 2:
        print("Error: Need at least two rows to determine the rule.")
        return
        
    width = len(grid[0])
    
    # 2. Deduce the rules from the grid evolution
    rules = {}
    
    # Iterate through grid rows (starting from the second) and columns
    for r in range(1, height):
        for c in range(width):
            # Determine the neighborhood from the previous row using periodic boundaries
            left_neighbor = grid[r-1][(c - 1 + width) % width]
            center_neighbor = grid[r-1][c]
            right_neighbor = grid[r-1][(c + 1) % width]
            
            neighborhood = (left_neighbor, center_neighbor, right_neighbor)
            
            # The result is the state of the current cell
            result = grid[r][c]
            
            # Store the rule if it's new
            if neighborhood not in rules:
                rules[neighborhood] = result

            # Once all 8 rules are found, we can stop
            if len(rules) == 8:
                break
        if len(rules) == 8:
            break

    # 3. Calculate the Wolfram code
    if len(rules) != 8:
        print(f"Error: Could not determine all 8 rules. Found {len(rules)}.")
        return

    # Standard order of neighborhoods for Wolfram code calculation
    neighborhood_order = [
        (1, 1, 1), (1, 1, 0), (1, 0, 1), (1, 0, 0),
        (0, 1, 1), (0, 1, 0), (0, 0, 1), (0, 0, 0)
    ]
    
    binary_string = "".join([str(rules[n]) for n in neighborhood_order])
    wolfram_code = int(binary_string, 2)

    # 4. Display the results
    print("Discovered Rule Set:")
    for neighborhood in neighborhood_order:
        print(f"  Neighborhood {neighborhood} -> Output {rules[neighborhood]}")
    
    print("\nConstructing Wolfram Code:")
    print(f"Binary representation (from rules for 111 down to 000): {binary_string}")
    
    calculation_parts = []
    for i, bit in enumerate(binary_string):
        power = 7 - i
        calculation_parts.append(f"{bit}*2^{power}")
    
    calculation_str = " + ".join(calculation_parts)
    print(f"Calculation: {calculation_str} = {wolfram_code}")
    
    print(f"\nThe Wolfram code for this cellular automaton is: {wolfram_code}")

solve_cellular_automaton()
<<<105>>>