import sys

def solve_cellular_automaton():
    """
    Identifies elementary cellular automaton rules consistent with the initial
    generations of the provided pattern.
    """
    # The grid represents the top 3 rows of the pattern (white=0, black=1).
    # We analyze only this part as later rows show inconsistencies.
    grid = [
        [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0]
    ]

    # The 8 neighborhoods in the standard order used for rule numbering.
    neighborhoods = [
        (1, 1, 1), (1, 1, 0), (1, 0, 1), (1, 0, 0),
        (0, 1, 1), (0, 1, 0), (0, 0, 1), (0, 0, 0)
    ]

    # This dictionary will store the observed output for each neighborhood.
    constraints = {}

    # Analyze the first two generations (rows 0->1 and 1->2).
    # We assume cells outside the grid are in state 0.
    for r in range(len(grid) - 1):
        row = grid[r]
        next_row = grid[r + 1]
        width = len(row)
        for c in range(width):
            left = row[c - 1] if c > 0 else 0
            center = row[c]
            right = row[c + 1] if c < width - 1 else 0
            
            hood = (left, center, right)
            output = next_row[c]
            
            # Store the constraint, ensuring no contradictions arise.
            if hood in constraints and constraints[hood] != output:
                print(f"Error: Contradictory rule found for neighborhood {hood}.", file=sys.stderr)
                return

            constraints[hood] = output

    # Build a template for the rule's 8-bit string, using 'x' for unconstrained parts.
    rule_template = []
    for hood in neighborhoods:
        if hood in constraints:
            rule_template.append(str(constraints[hood]))
        else:
            # This neighborhood was not observed, so its output is unconstrained.
            rule_template.append('x')
    
    rule_template_str = "".join(rule_template)

    # Find all valid rules by replacing 'x' with '0' and '1'.
    possible_rules = []
    def find_all_rules(template):
        if 'x' not in template:
            possible_rules.append(int(template, 2))
            return
        
        # Recurse for both possibilities of the first 'x'.
        find_all_rules(template.replace('x', '0', 1))
        find_all_rules(template.replace('x', '1', 1))

    find_all_rules(rule_template_str)
    
    # Sort the final list of rules.
    possible_rules.sort()
    
    # Print the result as a comma-separated string.
    print(','.join(map(str, possible_rules)))

solve_cellular_automaton()