import itertools

def solve_cellular_automaton():
    """
    Determines the possible elementary cellular automaton rules that could produce
    the behavior shown in the first three rows of the image.
    """
    # The 8 neighborhoods in the standard order used for rule numbering (111, 110, ..., 000).
    neighborhoods_ordered = list(itertools.product([1, 0], repeat=3))

    # Constraints are derived by observing the generation of rows 1 and 2.
    # R0: ...00100...
    # R1: ...01110...
    # R2: ...11011... (active part)

    # From R0 -> R1:
    # f(0,0,0) -> 0 (background)
    # f(0,0,1) -> 1 (left side of R1)
    # f(0,1,0) -> 1 (center of R1)
    # f(1,0,0) -> 1 (right side of R1)
    
    # From R1 -> R2:
    # f(0,1,1) -> 1 (left of center of R2)
    # f(1,1,1) -> 0 (center of R2)
    # f(1,1,0) -> 1 (right of center of R2)
    
    # The neighborhood (1,0,1) is not observed in generating the first three rows.
    # Its output (bit b5) is unconstrained.
    
    constraints = {
        (1,1,1): 0, # b7
        (1,1,0): 1, # b6
        # (1,0,1) -> ? (b5)
        (1,0,0): 1, # b4
        (0,1,1): 1, # b3
        (0,1,0): 1, # b2
        (0,0,1): 1, # b1
        (0,0,0): 0  # b0
    }

    possible_rules = []
    
    # Iterate through the two possibilities for the unconstrained neighborhood (1,0,1).
    for b5_value in [0, 1]:
        final_constraints = constraints.copy()
        final_constraints[(1,0,1)] = b5_value
        
        rule_binary_list = []
        for neighborhood in neighborhoods_ordered:
            rule_binary_list.append(str(final_constraints[neighborhood]))
            
        rule_binary_str = "".join(rule_binary_list)
        rule_decimal = int(rule_binary_str, 2)
        possible_rules.append(rule_decimal)

    possible_rules.sort()
    
    print("The possible rules are derived from the first three rows of the automaton's evolution.")
    print("The transitions observed are:")
    print("f(0,0,0) -> 0")
    print("f(0,0,1) -> 1")
    print("f(0,1,0) -> 1")
    print("f(1,0,0) -> 1")
    print("f(0,1,1) -> 1")
    print("f(1,1,1) -> 0")
    print("f(1,1,0) -> 1")
    print("\nThe neighborhood (1,0,1) is not observed, leading to two possibilities for its output.")
    print("\nCase 1: f(1,0,1) -> 0")
    print("Rule bits (b7..b0): 01011110")
    print("Decimal value: 64 + 16 + 8 + 4 + 2 = 94")
    
    print("\nCase 2: f(1,0,1) -> 1")
    print("Rule bits (b7..b0): 01111110")
    print("Decimal value: 64 + 32 + 16 + 8 + 4 + 2 = 126")

    print("\nFinal list of possible rules, sorted:")
    print(','.join(map(str, possible_rules)))

solve_cellular_automaton()