import collections

def solve_cellular_automaton():
    """
    Analyzes the provided cellular automaton pattern to determine its Wolfram code.
    """
    pattern_str = """
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
    # 1. Parse the pattern into a grid of 0s and 1s
    lines = pattern_str.strip().split('\n')
    grid = [[1 if char == '█' else 0 for char in line] for line in lines]
    height = len(grid)
    width = len(grid[0])

    # 2. Collect votes for each of the 8 possible rules
    rules_votes = collections.defaultdict(lambda: collections.defaultdict(int))

    for r in range(height - 1):
        for c in range(width):
            left_neighbor = grid[r][(c - 1 + width) % width]
            center_cell = grid[r][c]
            right_neighbor = grid[r][(c + 1) % width]
            neighborhood = (left_neighbor, center_cell, right_neighbor)
            outcome = grid[r + 1][c]
            rules_votes[neighborhood][outcome] += 1

    # 3. Determine the final rule set based on the votes
    ordered_neighborhoods = [
        (1, 1, 1), (1, 1, 0), (1, 0, 1), (1, 0, 0),
        (0, 1, 1), (0, 1, 0), (0, 0, 1), (0, 0, 0)
    ]
    binary_string = ""
    for neighborhood in ordered_neighborhoods:
        votes = rules_votes[neighborhood]
        if votes[1] > votes[0]:
            binary_string += "1"
        else:
            binary_string += "0"
            
    # 4. Calculate the Wolfram code from the binary string
    wolfram_code = int(binary_string, 2)

    # 5. Output the result, showing the calculation
    print(f"The best-fit binary rule string is: {binary_string}")
    print("\nThis binary string represents the outcomes for neighborhoods from 111 down to 000.")
    print("The Wolfram Code is its decimal equivalent.")

    number_parts = []
    for i, bit in enumerate(binary_string):
        power = 7 - i
        if bit == '1':
            value = 2**power
            number_parts.append(str(value))
    
    print("\nCalculation:")
    print(f"{' + '.join(number_parts)} = {wolfram_code}")

solve_cellular_automaton()
print("\n<<<45>>>")