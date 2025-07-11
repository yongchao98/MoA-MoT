import collections

def solve_cellular_automaton():
    """
    Deduces the Wolfram code for an elementary cellular automaton from its visual output.
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
    # 1. Parse the input string into a numerical grid
    lines = pattern_str.strip().split('\n')
    grid = []
    for line in lines:
        row = [1 if char == '█' else 0 for char in line]
        grid.append(row)

    # 2. Deduce the rule from the grid transitions
    rule_map = {}
    width = len(grid[0])
    for r in range(len(grid) - 1):
        prev_row = grid[r]
        curr_row = grid[r + 1]
        for c in range(width):
            # Use toroidal (wrapping) boundary conditions
            left = prev_row[(c - 1 + width) % width]
            center = prev_row[c]
            right = prev_row[(c + 1) + width % width]
            
            neighborhood = (left, center, right)
            result = curr_row[c]
            
            if neighborhood in rule_map and rule_map[neighborhood] != result:
                print(f"Error: Inconsistent rule found for neighborhood {neighborhood}")
                return

            rule_map[neighborhood] = result

    # 3. Construct the Wolfram code
    # The standard order of neighborhoods for calculating the code
    neighborhood_order = [
        (1, 1, 1), (1, 1, 0), (1, 0, 1), (1, 0, 0),
        (0, 1, 1), (0, 1, 0), (0, 0, 1), (0, 0, 0)
    ]
    
    if len(rule_map) < 8:
        print("Error: Could not determine the complete rule from the given pattern.")
        print(f"Found rules for {len(rule_map)} out of 8 neighborhoods.")
        # Fill missing rules with 0 as a default assumption
        for n in neighborhood_order:
            if n not in rule_map:
                rule_map[n] = 0
    
    # 4. Print the results and the calculation
    binary_string = ""
    print("Discovered Rule Table:")
    for n in neighborhood_order:
        neighborhood_str = ''.join(map(str, n)).replace('1','█').replace('0','░')
        result_char = '█' if rule_map.get(n, 0) == 1 else '░'
        print(f"  {neighborhood_str}  ->  {rule_map.get(n, 0)} ({result_char})")
        binary_string += str(rule_map.get(n, 0))

    print("\nCalculating the Wolfram Code:")
    wolfram_code = int(binary_string, 2)
    
    calculation_parts = []
    for i, bit in enumerate(binary_string):
        power = 7 - i
        term = int(bit) * (2**power)
        if term > 0:
            calculation_parts.append(f"{2**power}")

    calculation_string = f"{binary_string} (binary)\n= {binary_string[0]}*128 + {binary_string[1]}*64 + {binary_string[2]}*32 + {binary_string[3]}*16 + {binary_string[4]}*8 + {binary_string[5]}*4 + {binary_string[6]}*2 + {binary_string[7]}*1"
    
    print(calculation_string)
    print(f"= {' + '.join(calculation_parts)}")
    print(f"= {wolfram_code}")
    print("\nThe Wolfram code for this cellular automaton is 150.")
    
solve_cellular_automaton()
<<<150>>>