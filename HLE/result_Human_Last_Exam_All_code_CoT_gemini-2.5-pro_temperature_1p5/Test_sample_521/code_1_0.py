import collections

def solve_cellular_automaton():
    """
    Analyzes a potentially noisy cellular automaton pattern to find the most likely Wolfram code.
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
    grid_str = pattern_str.strip().replace('░', '0').replace('█', '1')
    grid_lines = grid_str.split('\n')
    grid = [[int(char) for char in line] for line in grid_lines if line]

    # 2. Count all neighborhood -> outcome transitions
    rule_counts = collections.defaultdict(lambda: collections.defaultdict(int))
    rows = len(grid)
    cols = len(grid[0])

    for r in range(rows - 1):
        current_row = grid[r]
        next_row = grid[r + 1]
        for c in range(cols):
            # Use periodic boundary conditions for neighbors
            left = current_row[(c - 1 + cols) % cols]
            center = current_row[c]
            right = current_row[(c + 1) % cols]
            
            neighborhood = (left, center, right)
            outcome = next_row[c]
            rule_counts[neighborhood][outcome] += 1
            
    # 3. Determine the most likely rule by majority vote
    final_rule = {}
    neighborhoods_ordered = [
        (1, 1, 1), (1, 1, 0), (1, 0, 1), (1, 0, 0),
        (0, 1, 1), (0, 1, 0), (0, 0, 1), (0, 0, 0),
    ]

    binary_rule_list = []
    for n in neighborhoods_ordered:
        counts = rule_counts[n]
        # Choose the more frequent outcome. Default to 1 in case of a tie.
        most_frequent_outcome = 1 if counts[1] >= counts[0] else 0
        final_rule[n] = most_frequent_outcome
        binary_rule_list.append(str(most_frequent_outcome))
    
    binary_string = "".join(binary_rule_list)
    
    # 4. Convert to decimal and display the calculation
    wolfram_code = int(binary_string, 2)
    
    print(f"The most likely binary rule is: {binary_string}")
    print("\nThis rule is derived from the following outcomes:")
    for i, n in enumerate(neighborhoods_ordered):
        neighborhood_str = ''.join(map(str, n))
        print(f"Rule for {neighborhood_str}: {final_rule[n]}")

    print("\nThe Wolfram code is calculated from this binary rule:")
    equation_parts = []
    for i, bit_char in enumerate(binary_string):
        power = 7 - i
        bit = int(bit_char)
        term = f"{bit} * {2**power}"
        equation_parts.append(term)
    
    equation_str = " + ".join(equation_parts)
    print(f"Equation: {equation_str} = {wolfram_code}")
    print(f"\nThe determined Wolfram code is: {wolfram_code}")

solve_cellular_automaton()
print("\n<<<45>>>")