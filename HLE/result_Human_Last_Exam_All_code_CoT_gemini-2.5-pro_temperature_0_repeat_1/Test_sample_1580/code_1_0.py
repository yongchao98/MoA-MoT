def solve_logic_puzzle():
    """
    Solves the Aric and Pi number puzzle by simulating their reasoning process.
    """
    aric_num = 13
    pi_num = 10
    possible_sums = [23, 20]

    # These sets will store the numbers that have been ruled out for each person.
    eliminated_for_aric = set()  # Numbers Aric can't have, from Pi's perspective
    eliminated_for_pi = set()    # Numbers Pi can't have, from Aric's perspective

    day = 1
    max_days = 20  # A safety break to prevent an infinite loop

    while day <= max_days:
        # --- Aric's Turn ---
        p_if_sum_23 = possible_sums[0] - aric_num
        p_if_sum_20 = possible_sums[1] - aric_num

        # Check if Aric can make a decision based on his current knowledge
        p1_is_possible = p_if_sum_23 not in eliminated_for_pi and p_if_sum_23 > 0
        p2_is_possible = p_if_sum_20 not in eliminated_for_pi and p_if_sum_20 > 0

        if p1_is_possible and not p2_is_possible:
            print(f"Aric announces the sum is {possible_sums[0]} on Day {day}.")
            print(f"The equation is: {aric_num} + {p_if_sum_23} = {possible_sums[0]}")
            return
        if not p1_is_possible and p2_is_possible:
            print(f"Aric announces the sum is {possible_sums[1]} on Day {day}.")
            print(f"The equation is: {aric_num} + {p_if_sum_20} = {possible_sums[1]}")
            return

        # Aric passes. Pi learns that Aric's number is not one that would have been decisive.
        # We calculate which hypothetical numbers for Aric would have led to a decision.
        decisive_a_numbers = set()
        for a_h in range(1, possible_sums[0]):
            if a_h in eliminated_for_aric:
                continue
            
            p1_h = possible_sums[0] - a_h
            p2_h = possible_sums[1] - a_h

            p1_h_possible = p1_h not in eliminated_for_pi and p1_h > 0
            p2_h_possible = p2_h not in eliminated_for_pi and p2_h > 0

            if p1_h_possible != p2_h_possible:  # XOR: one is possible but not the other
                decisive_a_numbers.add(a_h)
        eliminated_for_aric.update(decisive_a_numbers)

        # --- Pi's Turn ---
        a_if_sum_23 = possible_sums[0] - pi_num
        a_if_sum_20 = possible_sums[1] - pi_num

        # Check if Pi can make a decision
        a1_is_possible = a_if_sum_23 not in eliminated_for_aric and a_if_sum_23 > 0
        a2_is_possible = a_if_sum_20 not in eliminated_for_aric and a_if_sum_20 > 0

        if a1_is_possible and not a2_is_possible:
            print(f"Pi announces the sum is {possible_sums[0]} on Day {day}.")
            print(f"The equation is: {a_if_sum_23} + {pi_num} = {possible_sums[0]}")
            return
        if not a1_is_possible and a2_is_possible:
            print(f"Pi announces the sum is {possible_sums[1]} on Day {day}.")
            print(f"The equation is: {a_if_sum_20} + {pi_num} = {possible_sums[1]}")
            return

        # Pi passes. Aric learns that Pi's number is not one that would have been decisive.
        decisive_p_numbers = set()
        for p_h in range(1, possible_sums[0]):
            if p_h in eliminated_for_pi:
                continue

            a1_h = possible_sums[0] - p_h
            a2_h = possible_sums[1] - p_h

            a1_h_possible = a1_h not in eliminated_for_aric and a1_h > 0
            a2_h_possible = a2_h not in eliminated_for_aric and a2_h > 0

            if a1_h_possible != a2_h_possible:
                decisive_p_numbers.add(p_h)
        eliminated_for_pi.update(decisive_p_numbers)

        day += 1

    print("NEVER")

solve_logic_puzzle()