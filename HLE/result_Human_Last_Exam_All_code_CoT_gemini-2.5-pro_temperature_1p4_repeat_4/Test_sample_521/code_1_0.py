import collections

def solve_cellular_automaton():
    """
    Analyzes a visual representation of a cellular automaton to find its Wolfram code.
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
    width = len(grid[0])
    
    # 2. Count frequencies of all rule transitions
    # key: ((neighborhood_tuple), result_bit), value: count
    rule_counts = collections.defaultdict(int)

    for i in range(height - 1):
        for j in range(width):
            # Using periodic boundary conditions
            left_neighbor = grid[i][(j - 1 + width) % width]
            center = grid[i][j]
            right_neighbor = grid[i][(j + 1) + width) % width]
            
            neighborhood = (left_neighbor, center, right_neighbor)
            result = grid[i+1][j]
            
            rule_counts[(neighborhood, result)] += 1
            
    # 3. Determine the most likely rule by voting
    # Standard Wolfram order for neighborhoods
    neighborhoods_ordered = [
        (1, 1, 1), (1, 1, 0), (1, 0, 1), (1, 0, 0),
        (0, 1, 1), (0, 1, 0), (0, 0, 1), (0, 0, 0)
    ]
    
    binary_rule = []
    
    for n_tuple in neighborhoods_ordered:
        count_0 = rule_counts.get((n_tuple, 0), 0)
        count_1 = rule_counts.get((n_tuple, 1), 0)
        
        # Choose the outcome with the higher frequency
        result = 1 if count_1 > count_0 else 0
        binary_rule.append(str(result))
        
    binary_rule_str = "".join(binary_rule)
    
    # 4. Calculate Wolfram Code
    wolfram_code = int(binary_rule_str, 2)
    
    # 5. Output the results
    print("The 8 bits of the determined rule (for neighborhoods 111, 110, ..., 000) are:")
    # "output each number in the final equation" -> printing the bits
    print(' '.join(binary_rule))
    print(f"\nThis binary string '{binary_rule_str}' corresponds to the Wolfram code:")
    print(wolfram_code)
    
    # Return the final answer in the specified format
    print(f"\n<<<{wolfram_code}>>>")

solve_cellular_automaton()