def solve_cellular_automaton():
    """
    Finds the intermediate step in a 2-generation cellular automaton sequence
    by testing all 256 possible elementary rules.
    """
    top_row_str = "01101001"
    bottom_row_str = "10000111"

    top_row = [int(c) for c in top_row_str]
    bottom_row = [int(c) for c in bottom_row_str]
    n = len(top_row)

    def apply_rule(rule_outputs, row):
        """
        Applies a given rule to a row using periodic boundary conditions.
        The rule_outputs list maps neighborhoods 111, 110, ..., 000
        to the bits of the rule number.
        """
        new_row = [0] * n
        for i in range(n):
            # Get neighbors with periodic (wrapped) boundaries
            left = row[(i - 1 + n) % n]
            center = row[i]
            right = row[(i + 1) % n]

            # Convert the 3-cell neighborhood to an index from 7 (111) to 0 (000)
            neighborhood_index = 4 * left + 2 * center + 1 * right
            
            # The rule's output for this neighborhood.
            # Wolfram's convention: bit 0 of rule is for 111, bit 7 is for 000.
            new_row[i] = rule_outputs[7 - neighborhood_index]
        return new_row

    # Iterate through all 256 possible rules
    for rule_num in range(256):
        # Get the 8-bit binary representation of the rule
        rule_bin = format(rule_num, '08b')
        rule_outputs = [int(b) for b in rule_bin]

        # First generation: top row -> intermediate row
        intermediate_row = apply_rule(rule_outputs, top_row)

        # Second generation: intermediate row -> bottom row
        calculated_bottom_row = apply_rule(rule_outputs, intermediate_row)

        # Check if the result matches the given bottom row
        if calculated_bottom_row == bottom_row:
            intermediate_row_str = "".join(map(str, intermediate_row))
            print(top_row_str)
            print(intermediate_row_str)
            print(bottom_row_str)
            return

# Run the solver
solve_cellular_automaton()
<<<11011101>>>