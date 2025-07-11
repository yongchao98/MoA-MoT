def solve_logic_puzzle():
    """
    Solves the Aric and Pi number puzzle by simulating their reasoning process day by day.
    """
    aric_num = 13
    pi_num = 10
    possible_sums = [20, 23]

    # Sets of numbers for Aric or Pi that would have been solved on a previous turn
    known_unsolvable_for_A = set(range(1, 23))
    known_unsolvable_for_P = set(range(1, 23))
    
    # These numbers would be solved immediately by Aric on Day 1
    s_a_day1 = {a for a in range(1, 23) if not (20 - a > 0 and 23 - a > 0)}

    print(f"Aric's number: {aric_num}, Pi's number: {pi_num}")
    print(f"Possible sums: {possible_sums}\n")

    print(f"Thinking process:")
    print("-" * 20)

    # Main loop simulating the days
    for day in range(1, 10):
        # --- Aric's Turn ---
        print(f"Day {day}, Aric's turn:")
        
        # Aric checks if Pi passing on all previous turns eliminates one of his possibilities.
        # His possibilities for Pi's number are `p1` and `p2`.
        p1 = possible_sums[0] - aric_num
        p2 = possible_sums[1] - aric_num

        can_solve = False
        if p1 not in known_unsolvable_for_P:
            print(f"Aric knows Pi's number cannot be {p1}, because if it were, Pi would have already found a solution.")
            print(f"Therefore, Aric deduces that Pi's number must be {p2}.")
            sum_val = aric_num + p2
            can_solve = True
        elif p2 not in known_unsolvable_for_P:
            print(f"Aric knows Pi's number cannot be {p2}, because if it were, Pi would have already found a solution.")
            print(f"Therefore, Aric deduces that Pi's number must be {p1}.")
            sum_val = aric_num + p1
            can_solve = True
        
        if can_solve:
            print(f"\nConclusion: Aric gives the answer on Day {day}.")
            print(f"The equation is: {aric_num} + {pi_num} = {sum_val}")
            print(f"Final Answer Day: {day}")
            return

        print("Aric cannot decide. He passes.")
        
        # Update common knowledge based on Aric passing.
        # Find all numbers 'a' that would have been solved on this turn.
        newly_solved_by_A = {a for a in known_unsolvable_for_A if ((possible_sums[0] - a not in known_unsolvable_for_P) or (possible_sums[1] - a not in known_unsolvable_for_P))}
        known_unsolvable_for_A -= newly_solved_by_A
        

        # --- Pi's Turn ---
        print(f"Day {day}, Pi's turn:")
        
        # Pi checks if Aric passing on all previous turns (including this one) eliminates one of his possibilities.
        # His possibilities for Aric's number are `a1` and `a2`.
        a1 = possible_sums[0] - pi_num
        a2 = possible_sums[1] - pi_num
        
        can_solve = False
        if a1 not in known_unsolvable_for_A:
            print(f"Pi knows Aric's number cannot be {a1}, because if it were, Aric would have already found a solution.")
            print(f"Therefore, Pi deduces that Aric's number must be {a2}.")
            sum_val = pi_num + a2
            can_solve = True
        elif a2 not in known_unsolvable_for_A:
            print(f"Pi knows Aric's number cannot be {a2}, because if it were, Aric would have already found a solution.")
            print(f"Therefore, Pi deduces that Aric's number must be {a1}.")
            sum_val = pi_num + a1
            can_solve = True
            
        if can_solve:
            print(f"\nConclusion: Pi gives the answer on Day {day}.")
            print(f"The equation is: {aric_num} + {pi_num} = {sum_val}")
            print(f"Final Answer Day: {day}")
            return
            
        print("Pi cannot decide. He passes.\n")

        # Update common knowledge based on Pi passing.
        # Find all numbers 'p' that would have been solved on this turn.
        newly_solved_by_P = {p for p in known_unsolvable_for_P if ((possible_sums[0] - p not in known_unsolvable_for_A) or (possible_sums[1] - p not in known_unsolvable_for_A))}
        known_unsolvable_for_P -= newly_solved_by_P

solve_logic_puzzle()