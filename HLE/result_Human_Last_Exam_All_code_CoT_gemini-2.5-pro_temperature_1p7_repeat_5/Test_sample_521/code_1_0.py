import sys

def solve_automaton():
    """
    Determines the Wolfram code for a given elementary cellular automaton pattern
    by analyzing its transitions and using a majority vote to handle inconsistencies.
    """
    input_pattern = """
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
    lines = input_pattern.strip().split('\n')
    grid = [[1 if char == '█' else 0 for char in line] for line in lines]
    
    num_rows = len(grid)
    if num_rows == 0:
        print("Error: Empty grid provided.")
        return
        
    width = len(grid[0])
    
    neighborhood_patterns = [
        (1, 1, 1), (1, 1, 0), (1, 0, 1), (1, 0, 0),
        (0, 1, 1), (0, 1, 0), (0, 0, 1), (0, 0, 0)
    ]
    
    rule_counts = {p: {0: 0, 1: 0} for p in neighborhood_patterns}

    for t in range(num_rows - 1):
        for x in range(width):
            left = grid[t][(x - 1 + width) % width]
            center = grid[t][x]
            right = grid[t][(x + 1) % width]
            
            neighborhood = (left, center, right)
            outcome = grid[t+1][x]
            
            if neighborhood in rule_counts:
                rule_counts[neighborhood][outcome] += 1

    binary_string = ""
    for pattern in neighborhood_patterns:
        counts = rule_counts[pattern]
        if counts[1] > counts[0]:
            outcome_bit = 1
        else: # Default to 0 on tie or if 0 has more votes
            outcome_bit = 0
        binary_string += str(outcome_bit)
        
    wolfram_code = int(binary_string, 2)
    
    equation_parts = []
    for i, bit in enumerate(binary_string):
        power_of_2 = 2**(7-i)
        equation_parts.append(f"{bit}*{power_of_2}")
        
    final_equation = f"{' + '.join(equation_parts)} = {wolfram_code}"
    print(final_equation)


solve_automaton()
<<<137>>>