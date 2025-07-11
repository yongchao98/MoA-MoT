import itertools

def solve_automaton_puzzle():
    """
    Finds the unique intermediate step in a 2-step elementary cellular automaton evolution.
    """
    top_str = "01101001"
    bottom_str = "10000111"
    n = len(top_str)

    top = [int(c) for c in top_str]
    bottom = [int(c) for c in bottom_str]
    
    solutions = []

    # Iterate through all 2^8 possible binary strings for the middle row.
    for middle_candidate_tuple in itertools.product([0, 1], repeat=n):
        middle_candidate = list(middle_candidate_tuple)
        
        # A dictionary to store the deduced rule mapping.
        # Key: integer value of the 3-cell neighborhood (0-7)
        # Value: output bit (0 or 1)
        rule = {}
        is_consistent = True

        # Phase 1: Deduce rule from the top -> middle transition.
        for i in range(n):
            # Assume periodic boundary conditions.
            left = top[(i - 1 + n) % n]
            center = top[i]
            right = top[(i + 1) % n]
            
            # Convert neighborhood to its integer value (0-7).
            neighborhood_val = left * 4 + center * 2 + right * 1
            output = middle_candidate[i]

            # If this part of the rule is already defined, it must be consistent.
            if neighborhood_val in rule and rule[neighborhood_val] != output:
                is_consistent = False
                break
            rule[neighborhood_val] = output
        
        if not is_consistent:
            continue

        # Phase 2: Check for consistency and further deduce the rule from the middle -> bottom transition.
        for i in range(n):
            left = middle_candidate[(i - 1 + n) % n]
            center = middle_candidate[i]
            right = middle_candidate[(i + 1) % n]

            neighborhood_val = left * 4 + center * 2 + right * 1
            output = bottom[i]
            
            # If this part of the rule is already defined, it must match.
            if neighborhood_val in rule and rule[neighborhood_val] != output:
                is_consistent = False
                break
            rule[neighborhood_val] = output

        if is_consistent:
            # If we reach here, this middle_candidate is a valid solution.
            solutions.append("".join(map(str, middle_candidate)))

    # As per the problem, we expect a single unique solution.
    if len(solutions) == 1:
        print(solutions[0])
    else:
        # Fallback in case of unexpected results.
        print(f"Found {len(solutions)} solutions. The problem implies there should be only one.")
        for s in solutions:
            print(s)

solve_automaton_puzzle()