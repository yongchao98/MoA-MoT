def solve_puzzle():
    """
    Solves the Aric and Pi number puzzle by simulating their reasoning process.
    """
    aric_num = 13
    pi_num = 10
    possible_sums = [23, 20]

    # These bounds represent the common knowledge about the possible range of numbers.
    # Initially, we only know they are positive.
    known_bounds = {'a_min': 1, 'a_max': float('inf'), 'p_min': 1, 'p_max': float('inf')}

    day = 0
    while day < 10:  # Safety break for the loop
        day += 1
        
        # --- Aric's Turn ---
        print(f"--- Day {day}: Aric's Turn ---")
        
        # Aric's two scenarios for Pi's number
        p_scenario_1 = possible_sums[0] - aric_num  # Assumes sum is 23 -> P=10
        p_scenario_2 = possible_sums[1] - aric_num  # Assumes sum is 20 -> P=7

        # Check if knowledge eliminates a scenario for Aric
        p1_is_possible = known_bounds['p_min'] <= p_scenario_1 <= known_bounds['p_max']
        p2_is_possible = known_bounds['p_min'] <= p_scenario_2 <= known_bounds['p_max']

        if p1_is_possible and not p2_is_possible:
            print(f"Aric knows Pi's number cannot be {p_scenario_2} based on the accumulated knowledge.")
            print(f"He concludes Pi's number must be {p_scenario_1}, so the sum is {aric_num + p_scenario_1}.")
            print("Final equation:")
            print(f"{aric_num} + {p_scenario_1} = {aric_num + p_scenario_1}")
            print(f"\nThey give an answer on Day {day}.")
            return day
        elif not p1_is_possible and p2_is_possible:
            print(f"Aric knows Pi's number cannot be {p_scenario_1} based on the accumulated knowledge.")
            print(f"He concludes Pi's number must be {p_scenario_2}, so the sum is {aric_num + p_scenario_2}.")
            print("Final equation:")
            print(f"{aric_num} + {p_scenario_2} = {aric_num + p_scenario_2}")
            print(f"\nThey give an answer on Day {day}.")
            return day
        
        # If Aric cannot decide, he passes. We calculate what he reveals.
        # Aric would have known if his number `a` was in a range that made one of Pi's possibilities impossible given current knowledge.
        decisive_a_range = set()
        for a in range(known_bounds['a_min'], int(known_bounds['a_max']) + 1):
            p1 = possible_sums[0] - a
            p2 = possible_sums[1] - a
            p1_poss = known_bounds['p_min'] <= p1 <= known_bounds['p_max']
            p2_poss = known_bounds['p_min'] <= p2 <= known_bounds['p_max']
            if p1_poss != p2_poss:
                decisive_a_range.add(a)

        # Aric passing means his number is not in that set. We can update the bounds.
        if decisive_a_range:
            new_a_max = min(known_bounds['a_max'], min(decisive_a_range) - 1) if aric_num < min(decisive_a_range) else known_bounds['a_max']
            new_a_min = max(known_bounds['a_min'], max(decisive_a_range) + 1) if aric_num > max(decisive_a_range) else known_bounds['a_min']
            if new_a_max != known_bounds['a_max']:
                known_bounds['a_max'] = new_a_max
                print(f"Aric passes. This implies his number is not in {sorted(list(decisive_a_range))}. New common knowledge: Aric's number is at most {new_a_max}.")
            elif new_a_min != known_bounds['a_min']:
                known_bounds['a_min'] = new_a_min
                print(f"Aric passes. This implies his number is not in {sorted(list(decisive_a_range))}. New common knowledge: Aric's number is at least {new_a_min}.")
        else:
             print("Aric passes. This provides no new information on the number ranges.")


        # --- Pi's Turn ---
        print(f"--- Day {day}: Pi's Turn ---")
        
        # Pi's two scenarios for Aric's number
        a_scenario_1 = possible_sums[0] - pi_num  # Assumes sum is 23 -> A=13
        a_scenario_2 = possible_sums[1] - pi_num  # Assumes sum is 20 -> A=10
        
        a1_is_possible = known_bounds['a_min'] <= a_scenario_1 <= known_bounds['a_max']
        a2_is_possible = known_bounds['a_min'] <= a_scenario_2 <= known_bounds['a_max']
        
        if a1_is_possible and not a2_is_possible:
            print(f"Pi knows Aric's number cannot be {a_scenario_2} based on the accumulated knowledge.")
            print(f"She concludes Aric's number must be {a_scenario_1}, so the sum is {pi_num + a_scenario_1}.")
            print("Final equation:")
            print(f"{a_scenario_1} + {pi_num} = {pi_num + a_scenario_1}")
            print(f"\nThey give an answer on Day {day}.")
            return day
        elif not a1_is_possible and a2_is_possible:
            print(f"Pi knows Aric's number cannot be {a_scenario_1} based on the accumulated knowledge.")
            print(f"She concludes Aric's number must be {a_scenario_2}, so the sum is {pi_num + a_scenario_2}.")
            print("Final equation:")
            print(f"{a_scenario_2} + {pi_num} = {pi_num + a_scenario_2}")
            print(f"\nThey give an answer on Day {day}.")
            return day
        
        # If Pi cannot decide, she passes. We calculate what she reveals.
        decisive_p_range = set()
        for p in range(known_bounds['p_min'], int(known_bounds['p_max']) + 1 if known_bounds['p_max'] != float('inf') else 50):
            a1 = possible_sums[0] - p
            a2 = possible_sums[1] - p
            a1_poss = known_bounds['a_min'] <= a1 <= known_bounds['a_max']
            a2_poss = known_bounds['a_min'] <= a2 <= known_bounds['a_max']
            if a1_poss != a2_poss:
                decisive_p_range.add(p)
                
        if decisive_p_range:
            new_p_max = min(known_bounds['p_max'], min(decisive_p_range) - 1) if pi_num < min(decisive_p_range) else known_bounds['p_max']
            new_p_min = max(known_bounds['p_min'], max(decisive_p_range) + 1) if pi_num > max(decisive_p_range) else known_bounds['p_min']
            if new_p_max != known_bounds['p_max']:
                known_bounds['p_max'] = new_p_max
                print(f"Pi passes. This implies her number is not in {sorted(list(decisive_p_range))}. New common knowledge: Pi's number is at most {new_p_max}.")
            elif new_p_min != known_bounds['p_min']:
                known_bounds['p_min'] = new_p_min
                print(f"Pi passes. This implies her number is not in {sorted(list(decisive_p_range))}. New common knowledge: Pi's number is at least {new_p_min}.")
        else:
            print("Pi passes. This provides no new information on the number ranges.")

    print("They never give an answer.")
    return "NEVER"

# Execute the simulation
final_day = solve_puzzle()
print(f"\n<<<Day {final_day}>>>")
