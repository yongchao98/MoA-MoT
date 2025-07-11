import collections

def solve():
    """
    This function determines the possible elementary cellular automaton rules
    based on the initial, non-contradictory part of the evolution shown in the image.
    """
    
    # --- Step 1: Define constraints from the image ---
    # The image shows a pattern that, upon full inspection, contains contradictions
    # for any single elementary cellular automaton rule. For example, the neighborhood
    # (1, 1, 0) (black, black, white) yields a 1 (black) in the transition from
    # Gen 1 to Gen 2, but a 0 (white) from Gen 2 to Gen 3.
    #
    # We hypothesize that the intended question can be solved by considering only the
    # initial, unambiguous part of the evolution, where contradictions are absent.
    # The first three rows (Generations 0, 1, and 2) are clear and consistent.
    
    # We will determine the constraints imposed by the G0 -> G1 and G1 -> G2 transitions.
    # A value of `None` means the rule for that neighborhood is unconstrained
    # by this initial data.
    # Neighborhoods are tuples (left, center, right). White=0, Black=1.
    constraints = {
        # neighborhood: required_outcome
        (0, 0, 0): 0,    # Quiescent background stays white.
        (0, 0, 1): 1,    # Left edge of the triangle growth from G0->G1.
        (0, 1, 0): 1,    # Center of the triangle growth from G0->G1.
        (1, 0, 0): 1,    # Right edge of the triangle growth from G0->G1.
        (0, 1, 1): 1,    # Formation of G2's left-side structure from G1 (...B B . B B...).
        (1, 0, 1): None, # This neighborhood does not appear in G0 or G1.
        (1, 1, 0): 1,    # Formation of G2's right-side structure from G1.
        (1, 1, 1): 0     # Formation of G2's central white cell from G1's BBB.
    }

    # --- Step 2: Test all 256 rules against the constraints ---
    possible_rules = []
    for rule_number in range(256):
        is_consistent = True
        for neighborhood_tuple, required_outcome in constraints.items():
            if required_outcome is None:
                continue

            # The Wolfram rule number is a bitmask. The output for a neighborhood
            # (b2, b1, b0) is the bit at the position corresponding to the integer
            # value of that binary sequence.
            index = neighborhood_tuple[0] * 4 + neighborhood_tuple[1] * 2 + neighborhood_tuple[2]
            
            # Extract the bit for the current neighborhood from the rule number.
            rule_output = (rule_number >> index) & 1

            if rule_output != required_outcome:
                is_consistent = False
                break
        
        if is_consistent:
            possible_rules.append(rule_number)

    # --- Step 3: Print the result in the required format ---
    possible_rules.sort()
    
    # Let's derive the two numbers manually for explanation
    # The rule in binary (O7 O6 O5 O4 O3 O2 O1 O0) is constrained as follows:
    # O7=R(111)=0, O6=R(110)=1, O5=R(101)=?, O4=R(100)=1,
    # O3=R(011)=1, O2=R(010)=1, O1=R(001)=1, O0=R(000)=0.
    # Binary pattern: 01?11110
    
    # Case 1: ? = 0 --> 01011110_2 = 64+16+8+4+2 = 94
    # Case 2: ? = 1 --> 01111110_2 = 64+32+16+8+4+2 = 126
    # So the possible rules are 94 and 126.
    
    print(','.join(map(str, possible_rules)))

solve()