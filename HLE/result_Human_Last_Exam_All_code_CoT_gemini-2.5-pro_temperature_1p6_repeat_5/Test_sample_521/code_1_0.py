import sys

def solve_cellular_automaton():
    """
    Analyzes a textual representation of an elementary cellular automaton to determine its Wolfram code.
    """
    ca_visual = """
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
    # 1. Parse the input into a numerical grid
    lines = ca_visual.strip().split('\n')
    grid = [[1 if char == '█' else 0 for char in line] for line in lines]
    
    height = len(grid)
    if height < 2:
        print("Not enough generations to determine the rule.", file=sys.stderr)
        return
    width = len(grid[0])

    # 2. Infer the rules by observing transitions
    rules = {}
    
    for row_idx in range(height - 1):
        current_gen = grid[row_idx]
        next_gen = grid[row_idx + 1]
        for col_idx in range(width):
            # Using periodic boundary conditions
            left = current_gen[(col_idx - 1 + width) % width]
            center = current_gen[col_idx]
            right = current_gen[(col_idx + 1) % width]
            
            neighborhood = (left, center, right)
            output = next_gen[col_idx]
            
            if neighborhood in rules and rules[neighborhood] != output:
                print(f"Inconsistent rule found for neighborhood {neighborhood}", file=sys.stderr)
                return
            rules[neighborhood] = output

    # 3. Construct the Wolfram code
    neighborhood_order = [
        (1, 1, 1), (1, 1, 0), (1, 0, 1), (1, 0, 0),
        (0, 1, 1), (0, 1, 0), (0, 0, 1), (0, 0, 0)
    ]

    if len(rules) < 8:
        print("Could not determine all 8 rules from the given pattern.", file=sys.stderr)
        return

    binary_string = ""
    for n in neighborhood_order:
        binary_string += str(rules[n])
        print(f"The rule for the neighborhood ({n[0]}, {n[1]}, {n[2]}) is {rules[n]}.")

    # 4. Calculate and display the final result
    wolfram_code = int(binary_string, 2)
    
    equation_parts = []
    powers_of_2 = [128, 64, 32, 16, 8, 4, 2, 1]
    bits = [int(b) for b in binary_string]
    
    for i in range(8):
        equation_parts.append(f"{powers_of_2[i]} * {bits[i]}")
        
    equation_str = " + ".join(equation_parts)
    
    print("\nThe binary representation of the rule is " + " ".join(list(binary_string)) + ".")
    print(f"The Wolfram code is calculated as: {equation_str} = {wolfram_code}")

solve_cellular_automaton()
<<<75>>>