import itertools

def solve():
    """
    Solves for the possible elementary cellular automaton rules based on observed transitions.
    """

    # Manually transcribed first 3 rows of the pattern, which are the clearest.
    # 0 = white, 1 = black. Width is cropped for simplicity.
    pattern = [
        [0, 0, 0, 1, 0, 0, 0],  # t=0
        [0, 0, 1, 1, 1, 0, 0],  # t=1
        [0, 1, 1, 0, 1, 1, 0],  # t=2
    ]

    # Deduce constraints from the pattern
    # A constraint is a mapping from a neighborhood (3-tuple) to an output (0 or 1)
    constraints = {}
    # Loop through the first two transitions (t=0 -> t=1 and t=1 -> t=2)
    for t in range(len(pattern) - 1):
        parent_row = pattern[t]
        child_row = pattern[t+1]
        # Iterate over cells where the neighborhood is fully contained within the row
        for i in range(1, len(parent_row) - 1):
            neighborhood = tuple(parent_row[i-1:i+2])
            output = child_row[i]
            
            # If we see a neighborhood for the first time, add it as a constraint
            if neighborhood not in constraints:
                constraints[neighborhood] = output
            # If we see it again, it must be consistent with the previous observation
            elif constraints[neighborhood] != output:
                print(f"Error: Inconsistent pattern found for neighborhood {neighborhood}.")
                return

    print("Deduced constraints from the first two transitions:")
    # Sort constraints for deterministic printing
    sorted_neighborhoods = sorted(constraints.keys(), reverse=True)
    for neighborhood in sorted_neighborhoods:
      print(f"Neighborhood {neighborhood} -> Output {constraints[neighborhood]}")
    
    # Identify which neighborhood rule is unconstrained
    all_neighborhoods = [tuple(p) for p in itertools.product([1, 0], repeat=3)]
    unconstrained = []
    for n in all_neighborhoods:
      if n not in constraints:
        unconstrained.append(n)
    
    if unconstrained:
      print(f"\nThe rule for neighborhood {unconstrained[0]} is unconstrained by these observations.")
    else:
      print("\nAll neighborhood rules are constrained.")


    # Find all rules (0-255) that satisfy the constraints
    matching_rules = []
    for rule_number in range(256):
        # The 8 bits of the rule number correspond to outputs for neighborhoods 111, 110, ..., 000
        rule_bits = format(rule_number, '08b')
        is_match = True
        for i, neighborhood_tuple in enumerate(all_neighborhoods):
            # Check against deduced constraints
            if neighborhood_tuple in constraints:
                # The index `i` maps neighborhood value (e.g., '111'=7) to rule bit position.
                # `rule_bits` is ordered from '111' (index 0) to '000' (index 7).
                # `all_neighborhoods` is created in ascending value (000..111)
                # so we use 7-i to get the correct index for wolfram code order
                rule_output = int(rule_bits[7-i])
                if constraints[neighborhood_tuple] != rule_output:
                    is_match = False
                    break
        
        if is_match:
            matching_rules.append(rule_number)

    print("\nThe following rules are consistent with the observed pattern:")
    print(",".join(map(str, sorted(matching_rules))))

solve()