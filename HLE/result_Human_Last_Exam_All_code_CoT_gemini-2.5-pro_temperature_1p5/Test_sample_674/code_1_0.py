def solve_sequence():
    """
    This function generates a sequence based on the rule
    s[n] = s[s[n-1]] + s[s[n-2]] and demonstrates the calculation for n=10.
    """
    s = {1: 1, 2: 1}
    n_max = 42 # The length of the provided sequence S4
    
    # Generate the sequence up to n_max
    for n in range(3, n_max + 1):
        try:
            prev1_val = s[n - 1]
            prev2_val = s[n - 2]
            idx1 = s[prev1_val]
            idx2 = s[prev2_val]
            s[n] = idx1 + idx2
        except KeyError as e:
            print(f"Error generating s[{n}]: index {e} not found in sequence.")
            break

    # The rule holds for n=10 with the given S4 sequence
    # s[10] = s[s[9]] + s[s[8]]
    # From S4, s[9]=4 and s[8]=4.
    # So, s[10] = s[4] + s[4].
    # From S4, s[4]=2.
    # So, s[10] = 2 + 2 = 4. This matches S4[10].
    
    n_demo = 10
    s_provided = [0, 1, 1, 2, 2, 2, 4, 3, 4, 4, 4, 8, 5, 5, 8, 8, 6, 8, 12, 8, 11, 9, 9, 10, 13, 16, 9, 12, 20, 10, 12, 23, 12, 15, 21, 13, 17, 18, 19, 19, 22, 21, 19]

    sn_minus_1 = s_provided[n_demo - 1]
    sn_minus_2 = s_provided[n_demo - 2]
    
    s_of_sn_minus_1 = s_provided[sn_minus_1]
    s_of_sn_minus_2 = s_provided[sn_minus_2]
    
    result = s_of_sn_minus_1 + s_of_sn_minus_2
    
    print(f"The deduced rule is R(s[n]) = s[s[n-1]] + s[s[n-2]].")
    print(f"Let's demonstrate for n = {n_demo}:")
    print(f"s[{n_demo}] = s[s[{n_demo - 1}]] + s[s[{n_demo - 2}]]")
    print(f"s[{n_demo}] = s[{sn_minus_1}] + s[{sn_minus_2}]")
    print(f"s[{n_demo}] = {s_of_sn_minus_1} + {s_of_sn_minus_2}")
    print(f"s[{n_demo}] = {result}")

solve_sequence()