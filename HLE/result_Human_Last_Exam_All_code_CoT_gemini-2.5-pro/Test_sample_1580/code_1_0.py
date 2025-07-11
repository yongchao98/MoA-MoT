def solve_logic_puzzle():
    """
    Solves the Aric and Pi number puzzle by simulating their reasoning process.
    """
    aric_num = 13
    pi_num = 10
    possible_sums = [20, 23]

    # Sets of numbers ruled out for each person based on their passes.
    excluded_for_aric = set()
    excluded_for_pi = set()

    # Max days to prevent infinite loop, though it should resolve quickly.
    for day in range(1, 20):
        # --- Aric's Turn ---
        
        # Aric determines his two possibilities for Pi's number.
        p_if_sum_is_20 = possible_sums[0] - aric_num
        p_if_sum_is_23 = possible_sums[1] - aric_num

        # He checks if the accumulated knowledge rules out one possibility.
        # A possibility is ruled out if the number is not positive or if Pi's
        # previous passes have shown that Pi cannot have that number.
        possibility_20_is_invalid = (p_if_sum_is_20 <= 0 or p_if_sum_is_20 in excluded_for_pi)
        possibility_23_is_invalid = (p_if_sum_is_23 <= 0 or p_if_sum_is_23 in excluded_for_pi)

        # If exactly one possibility is invalid, Aric knows the answer.
        if possibility_20_is_invalid and not possibility_23_is_invalid:
            print(f"Aric answers on Day {day}.")
            print(f"{aric_num} + {p_if_sum_is_23} = {aric_num + p_if_sum_is_23}")
            return
        if possibility_23_is_invalid and not possibility_20_is_invalid:
            print(f"Aric answers on Day {day}.")
            print(f"{aric_num} + {p_if_sum_is_20} = {aric_num + p_if_sum_is_20}")
            return

        # If Aric passes, we identify which numbers he would have had to solve the puzzle.
        # Pi will learn that Aric does not have one of these numbers.
        newly_excluded_for_aric = set()
        for a_candidate in range(1, possible_sums[1]):
            if a_candidate in excluded_for_aric:
                continue
            
            p_cand_1 = possible_sums[0] - a_candidate
            p_cand_2 = possible_sums[1] - a_candidate

            cand1_impossible = p_cand_1 <= 0 or p_cand_1 in excluded_for_pi
            cand2_impossible = p_cand_2 <= 0 or p_cand_2 in excluded_for_pi

            if cand1_impossible ^ cand2_impossible: # XOR: one is impossible, but not both
                newly_excluded_for_aric.add(a_candidate)
        excluded_for_aric.update(newly_excluded_for_aric)

        # --- Pi's Turn ---
        
        # Pi determines his two possibilities for Aric's number.
        a_if_sum_is_20 = possible_sums[0] - pi_num
        a_if_sum_is_23 = possible_sums[1] - pi_num

        # He checks if Aric's passes have ruled out one possibility.
        possibility_20_is_invalid = (a_if_sum_is_20 <= 0 or a_if_sum_is_20 in excluded_for_aric)
        possibility_23_is_invalid = (a_if_sum_is_23 <= 0 or a_if_sum_is_23 in excluded_for_aric)

        # If exactly one possibility is invalid, Pi knows the answer.
        if possibility_20_is_invalid and not possibility_23_is_invalid:
            print(f"Pi answers on Day {day}.")
            # The problem asks to output the numbers in the equation
            print(f"{a_if_sum_is_23} + {pi_num} = {a_if_sum_is_23 + pi_num}")
            return
        if possibility_23_is_invalid and not possibility_20_is_invalid:
            print(f"Pi answers on Day {day}.")
            print(f"{a_if_sum_is_20} + {pi_num} = {a_if_sum_is_20 + pi_num}")
            return

        # If Pi passes, we update the set of numbers that Aric will know Pi cannot have.
        newly_excluded_for_pi = set()
        for p_candidate in range(1, possible_sums[1]):
            if p_candidate in excluded_for_pi:
                continue
            
            a_cand_1 = possible_sums[0] - p_candidate
            a_cand_2 = possible_sums[1] - p_candidate

            cand1_impossible = a_cand_1 <= 0 or a_cand_1 in excluded_for_aric
            cand2_impossible = a_cand_2 <= 0 or a_cand_2 in excluded_for_aric

            if cand1_impossible ^ cand2_impossible:
                newly_excluded_for_pi.add(p_candidate)
        excluded_for_pi.update(newly_excluded_for_pi)

    # If the loop completes, they never find a certain answer.
    print("NEVER")

if __name__ == '__main__':
    solve_logic_puzzle()