import itertools

def find_ca_rules():
    """
    Analyzes a cellular automaton pattern to determine the generating rule(s).
    """
    # Step 1: Represent the grid from the image (White=0, Black=1).
    grid = [
        [0,0,0,0,0,0,0,1,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,1,1,1,0,0,0,0,0,0],
        [0,0,0,0,0,1,1,0,1,1,0,0,0,0,0],
        [0,0,0,0,1,1,0,0,0,1,1,0,0,0,0],
        [0,0,0,1,1,0,1,0,1,0,1,1,0,0,0],
        [0,0,1,1,0,1,1,1,1,1,0,1,1,0,0],
        [0,1,1,0,1,1,0,0,0,1,1,0,1,1,0],
        [1,1,0,1,1,0,1,0,1,0,1,1,0,1,1]
    ]

    # Step 2 & 3: Extract transitions and check for contradictions.
    # We use a heuristic: the first observed transition for a neighborhood is assumed correct.
    constraints = {}
    contradictions_found = False
    
    num_rows = len(grid)
    num_cols = len(grid[0])

    for t in range(num_rows - 1):
        for j in range(1, num_cols - 1):
            neighborhood = tuple(grid[t][j-1 : j+2])
            child = grid[t+1][j]
            
            if neighborhood not in constraints:
                constraints[neighborhood] = child
            elif constraints[neighborhood] != child:
                # A contradiction is found, but we proceed with the heuristic.
                contradictions_found = True

    if contradictions_found:
        print("Warning: Inconsistencies were found in the image pattern.")
        print("The same neighborhood pattern produces different outcomes.")
        print("Proceeding with the heuristic that the first observed outcome for each neighborhood is correct.")
        print("-" * 20)

    # Step 4 & 5: Identify unobserved neighborhoods and generate all possible rules.
    all_neighborhoods = [tuple(i) for i in itertools.product([1, 0], repeat=3)]
    unobserved = [n for n in all_neighborhoods if n not in constraints]

    possible_rules = []
    
    # Generate all combinations for the "don't care" outputs.
    # If there are no unobserved neighborhoods, this loop runs once.
    for combo in itertools.product([0, 1], repeat=len(unobserved)):
        # Create a complete rule set for this combination
        final_rule = constraints.copy()
        for i, n in enumerate(unobserved):
            final_rule[n] = combo[i]
        
        # Build the 8-bit binary string for the rule number
        rule_binary_str = ""
        for n in all_neighborhoods:
            rule_binary_str += str(final_rule[n])
            
        rule_decimal = int(rule_binary_str, 2)
        possible_rules.append(rule_decimal)

    # Step 6: Sort and output the result.
    possible_rules.sort()
    
    print("Based on the analysis, the following rule(s) could have produced the behavior (under the stated heuristic if inconsistencies were found):")
    # The final print is just the numbers as requested.
    print(','.join(map(str, possible_rules)))

find_ca_rules()