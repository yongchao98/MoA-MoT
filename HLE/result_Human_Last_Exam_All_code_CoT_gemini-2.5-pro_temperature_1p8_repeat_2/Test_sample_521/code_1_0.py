import sys

def solve_cellular_automaton():
    """
    Analyzes a given cellular automaton pattern to determine its Wolfram code.
    """
    
    # The visual pattern of the cellular automaton's evolution.
    pattern_string = """
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
    # Convert pattern to binary representation (1 for '█', 0 for '░')
    lines = [line.strip().replace('█', '1').replace('░', '0') for line in pattern_string.strip().split('\n')]
    if not lines or not all(len(line) == len(lines[0]) for line in lines):
        print("Error: The pattern is empty or has lines of inconsistent lengths.")
        return

    width = len(lines[0])
    rules = {}

    # Iterate through consecutive lines to deduce the rules
    for t in range(len(lines) - 1):
        if len(rules) == 8:
            break  # Stop once all 8 rules are found
        
        current_line = lines[t]
        next_line = lines[t + 1]
        
        # Iterate through cells, assuming periodic boundary conditions
        for i in range(width):
            left = int(current_line[(i - 1 + width) % width])
            center = int(current_line[i])
            right = int(current_line[(i + 1) % width])
            
            neighborhood = (left, center, right)
            result = int(next_line[i])

            if neighborhood not in rules:
                rules[neighborhood] = result
            elif rules[neighborhood] != result:
                # This would indicate an inconsistency, but we assume the pattern is valid.
                pass

    # The standard order for neighborhood evaluation
    neighborhood_order = [
        (1, 1, 1), (1, 1, 0), (1, 0, 1), (1, 0, 0),
        (0, 1, 1), (0, 1, 0), (0, 0, 1), (0, 0, 0)
    ]

    if len(rules) != 8:
        print("Could not deduce all 8 rules from the provided pattern.")
        return

    binary_rule_str = ""
    for neighborhood in neighborhood_order:
        binary_rule_str += str(rules[neighborhood])
    
    wolfram_code = int(binary_rule_str, 2)
    
    # Construct the calculation string
    calculation_parts = []
    total = 0
    for i, bit in enumerate(binary_rule_str):
        power_of_2 = 2**(7 - i)
        term_value = int(bit) * power_of_2
        total += term_value
        calculation_parts.append(f"{bit}*{power_of_2}")
    
    calculation_string = " + ".join(calculation_parts)

    print(f"The automaton follows Rule {wolfram_code}.")
    print("\nThe binary representation of the rule is:", binary_rule_str)
    print("\nThis is calculated from the outcomes of the 8 neighborhoods:")
    for neighborhood, result in zip(neighborhood_order, binary_rule_str):
        print(f"  {neighborhood[0]}{neighborhood[1]}{neighborhood[2]} -> {result}")

    print(f"\nThe Wolfram code is the decimal value of this binary number:")
    print(f"{calculation_string} = {total}")

solve_cellular_automaton()
<<<77>>>