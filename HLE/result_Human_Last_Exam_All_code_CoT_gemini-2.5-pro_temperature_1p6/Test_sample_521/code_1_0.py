import textwrap

def solve_wolfram_code():
    """
    Analyzes the provided cellular automaton evolution to determine its Wolfram code.
    """
    # The visual representation of the cellular automaton's evolution.
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

    # 1. Parse the visual input into a numerical grid (0s and 1s).
    lines = textwrap.dedent(ca_visual).strip().split('\n')
    grid = []
    for line in lines:
        row = [1 if char == '█' else 0 for char in line.strip()]
        grid.append(row)

    # 2. Deduce the rules by observing transitions.
    rules = {}
    rows, cols = len(grid), len(grid[0])
    for r in range(rows - 1):
        for c in range(1, cols - 1):
            neighborhood = tuple(grid[r][c-1:c+2])
            result = grid[r+1][c]
            if neighborhood not in rules:
                rules[neighborhood] = result
            # Stop once all 8 rules are found.
            if len(rules) == 8:
                break
        if len(rules) == 8:
            break

    # 3. Construct the Wolfram code from the discovered rules.
    # The standard order for neighborhoods in Wolfram code calculation.
    neighborhood_order = [
        (1, 1, 1), (1, 1, 0), (1, 0, 1), (1, 0, 0),
        (0, 1, 1), (0, 1, 0), (0, 0, 1), (0, 0, 0)
    ]

    binary_string = ""
    for n in neighborhood_order:
        binary_string += str(rules[n])

    wolfram_code = int(binary_string, 2)

    # 4. Print the detailed calculation and the final result.
    print("Discovered Rules:")
    for n in neighborhood_order:
        print(f"Rule for {''.join(map(str, n))} is {rules[n]}")
    
    print("\nBinary string (from 111 down to 000):")
    print(binary_string)

    print("\nCalculating Wolfram Code:")
    calculation_parts = []
    for i, bit in enumerate(binary_string):
        power_of_2 = 2**(7 - i)
        calculation_parts.append(f"{bit}*{power_of_2}")
    
    calculation_str = " + ".join(calculation_parts)
    print(f"Wolfram code = {calculation_str} = {wolfram_code}")

    print(f"\nThe Wolfram code of this cellular automaton is:")
    print(wolfram_code)


solve_wolfram_code()
<<<41>>>