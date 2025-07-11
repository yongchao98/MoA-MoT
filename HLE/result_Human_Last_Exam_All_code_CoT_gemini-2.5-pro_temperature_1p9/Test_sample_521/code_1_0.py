import sys

def find_wolfram_code():
    """
    Analyzes the evolution of an elementary cellular automaton to determine its
    Wolfram code.
    """
    ca_evolution = """
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
    # 1. Parse into a grid of 0s and 1s
    lines = ca_evolution.strip().split('\n')
    grid = [[1 if char == '█' else 0 for char in line] for line in lines if line]

    if len(grid) < 2 or not grid[0]:
        print("Error: Input data is insufficient to determine the rule.", file=sys.stderr)
        return

    num_rows, width = len(grid), len(grid[0])
    rules = {}

    # 2. Extract rules by observing transitions between rows
    for r in range(num_rows - 1):
        current_row = grid[r]
        next_row = grid[r+1]
        for c in range(width):
            # Determine neighborhood using periodic (wraparound) boundary conditions
            left = current_row[(c - 1 + width) % width]
            center = current_row[c]
            right = current_row[(c + 1) % width]
            
            neighborhood = (left, center, right)
            result = next_row[c]

            # Store the observed rule. Check for contradictions.
            if neighborhood in rules and rules[neighborhood] != result:
                print(f"Error: Contradiction found for rule {neighborhood}.", file=sys.stderr)
                return
            rules[neighborhood] = result

    # 3. Reconstruct the 8-bit binary string for the Wolfram code
    neighborhood_order = [
        (1, 1, 1), (1, 1, 0), (1, 0, 1), (1, 0, 0),
        (0, 1, 1), (0, 1, 0), (0, 0, 1), (0, 0, 0)
    ]
    
    binary_rule_list = []
    for neighborhood in neighborhood_order:
        if neighborhood not in rules:
            print(f"Error: Rule for neighborhood {neighborhood} not found.", file=sys.stderr)
            return
        binary_rule_list.append(str(rules[neighborhood]))
        
    binary_rule_string = "".join(binary_rule_list)
    
    # 4. Convert the binary string to its decimal equivalent
    wolfram_code = int(binary_rule_string, 2)
    
    # 5. Output the result, showing the numbers that form the final calculation
    print("The rule set was inferred from the provided grid.")
    print(f"The 8-bit binary representation of the rule is: {binary_rule_string}")
    print(f"The decimal value of {binary_rule_string} is {wolfram_code}.")
    print("\nTherefore, the Wolfram code of this cellular automaton is:")
    print(wolfram_code)

if __name__ == "__main__":
    find_wolfram_code()