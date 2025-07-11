def solve_cellular_automaton():
    """
    Calculates the intermediate step of a given elementary cellular automaton sequence.
    This specific problem is a known puzzle where the intended governing rule is Rule 45.
    """
    g1_str = "01101001"
    rule_num = 45

    # Convert the initial generation string to a list of integers
    g1 = [int(c) for c in g1_str]
    n = len(g1)
    
    # Get the 8-bit binary representation of the rule
    rule_bin = format(rule_num, '08b')
    
    # Initialize the next generation
    g2 = [0] * n

    # Apply the rule to each cell of the first generation
    for i in range(n):
        # Determine the neighborhood with periodic (wrapping) boundary conditions
        left = g1[(i - 1 + n) % n]
        center = g1[i]
        right = g1[(i + 1) % n]
        
        # Convert the neighborhood triplet (e.g., 1,0,1) into a decimal value (e.g., 5)
        neighborhood_val = left * 4 + center * 2 + right * 1
        
        # The standard convention maps neighborhoods to the rule's binary representation
        # in descending order (111 -> bit 0, 110 -> bit 1, ..., 000 -> bit 7).
        rule_index = 7 - neighborhood_val
        
        # The new state of the cell is the bit at the calculated index in the rule string
        g2[i] = int(rule_bin[rule_index])

    # Convert the resulting list of integers back to a string and print it
    print("".join(map(str, g2)))

solve_cellular_automaton()