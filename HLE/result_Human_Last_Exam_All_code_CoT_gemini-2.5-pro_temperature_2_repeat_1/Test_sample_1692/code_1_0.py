def solve_cellular_automaton():
    """
    Finds the intermediate step in a 2-step cellular automaton evolution by testing all 256 possible rules.
    """
    top_row_str = "01101001"
    bottom_row_str = "10000111"

    top_row = [int(c) for c in top_row_str]
    bottom_row = [int(c) for c in bottom_row_str]
    n = len(top_row)

    # Iterate through all 256 possible elementary rules
    for rule_number in range(256):
        # The rule's binary representation defines the output for each of the 8 possible neighbor triplets
        # The standard Wolfram convention orders triplets from '111' down to '000'
        rule_bin = format(rule_number, '08b')
        rule_map = {
            (1, 1, 1): int(rule_bin[0]),
            (1, 1, 0): int(rule_bin[1]),
            (1, 0, 1): int(rule_bin[2]),
            (1, 0, 0): int(rule_bin[3]),
            (0, 1, 1): int(rule_bin[4]),
            (0, 1, 0): int(rule_bin[5]),
            (0, 0, 1): int(rule_bin[6]),
            (0, 0, 0): int(rule_bin[7]),
        }

        # Step 1: Calculate the intermediate row from the top row
        intermediate_row = [0] * n
        for i in range(n):
            # Get neighbors using periodic boundary conditions
            left = top_row[(i - 1) % n]
            center = top_row[i]
            right = top_row[(i + 1) % n]
            triplet = (left, center, right)
            intermediate_row[i] = rule_map[triplet]

        # Step 2: Calculate the bottom row from the intermediate row
        calculated_bottom_row = [0] * n
        for i in range(n):
            # Get neighbors using periodic boundary conditions
            left = intermediate_row[(i - 1) % n]
            center = intermediate_row[i]
            right = intermediate_row[(i + 1) % n]
            triplet = (left, center, right)
            calculated_bottom_row[i] = rule_map[triplet]
        
        # Check if the calculated bottom row matches the target bottom row
        if calculated_bottom_row == bottom_row:
            intermediate_row_str = "".join(map(str, intermediate_row))
            print(top_row_str)
            print(intermediate_row_str)
            print(bottom_row_str)
            # This is the unique solution
            print(f"<<<{intermediate_row_str}>>>")
            return

solve_cellular_automaton()