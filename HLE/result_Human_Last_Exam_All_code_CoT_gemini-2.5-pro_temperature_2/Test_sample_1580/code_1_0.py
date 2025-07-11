import math

def solve_puzzle():
    """
    Simulates the logic puzzle to determine which day the sum is known.
    """
    ARIC_NUM = 13
    PI_NUM = 10
    SUMS = [20, 23]
    MAX_DAYS = 20

    # Sets of numbers that would have been revealed on previous turns
    known_p_numbers = set()
    known_a_numbers = set()

    # Public information gained from passes
    # Aric's first pass reveals a < 20. Pi's pass reveals p < 20.
    # We build on top of this by adding numbers revealed through more complex logic.
    initial_a_knowledge = {a for a in range(SUMS[0], SUMS[1])}
    initial_p_knowledge = {p for p in range(SUMS[0], SUMS[1])}
    
    known_a_numbers.update(initial_a_knowledge)
    known_p_numbers.update(initial_p_knowledge)

    for day in range(1, MAX_DAYS + 1):
        # --- Aric's Turn (Morning) ---
        
        # Aric checks if Pi passing on all previous turns eliminates one of his possibilities.
        p_possibilities = {SUMS[0] - ARIC_NUM, SUMS[1] - ARIC_NUM}
        unresolved_p = p_possibilities - known_p_numbers
        
        if len(unresolved_p) == 1:
            winning_p = unresolved_p.pop()
            eliminated_p = p_possibilities.difference({winning_p}).pop()
            
            print(f"On Day {day}, Aric knows the answer.")
            print(f"He considered two possibilities for Pi's number: {p_possibilities}.")
            print(f"From previous rounds of passing, he knew Pi's number could not be in the set {sorted(list(known_p_numbers))}.")
            print(f"This information eliminates {eliminated_p} as a possibility for Pi's number.")
            print(f"Therefore, Pi's number must be {winning_p}.")
            final_sum = ARIC_NUM + winning_p
            print(f"The final equation is: {ARIC_NUM} + {winning_p} = {final_sum}")
            return

        # If Aric passes, we update the public knowledge of what Aric's number *cannot* be.
        # These are the numbers that *would* have solved it this turn.
        newly_solved_a = set()
        for a_test in range(1, SUMS[1]):
            if a_test in known_a_numbers:
                continue
            a_test_p_poss = {SUMS[0] - a_test, SUMS[1] - a_test}
            if len(a_test_p_poss.intersection(known_p_numbers)) == 1:
                newly_solved_a.add(a_test)
        known_a_numbers.update(newly_solved_a)

        # --- Pi's Turn (Afternoon) ---

        # Pi checks if Aric passing on this and all previous turns eliminates one of her possibilities.
        a_possibilities = {SUMS[0] - PI_NUM, SUMS[1] - PI_NUM}
        unresolved_a = a_possibilities - known_a_numbers
        
        if len(unresolved_a) == 1:
            winning_a = unresolved_a.pop()
            eliminated_a = a_possibilities.difference({winning_a}).pop()

            print(f"On Day {day}, Pi knows the answer.")
            print(f"She considered two possibilities for Aric's number: {a_possibilities}.")
            print(f"From previous rounds of passing, she knew Aric's number could not be in the set {sorted(list(known_a_numbers))}.")
            print(f"This information eliminates {eliminated_a} as a possibility for Aric's number.")
            print(f"Therefore, Aric's number must be {winning_a}.")
            final_sum = PI_NUM + winning_a
            print(f"The final equation is: {winning_a} + {PI_NUM} = {final_sum}")
            return
            
        # If Pi passes, update public knowledge about Pi's number
        newly_solved_p = set()
        for p_test in range(1, SUMS[1]):
            if p_test in known_p_numbers:
                continue
            p_test_a_poss = {SUMS[0] - p_test, SUMS[1] - p_test}
            if len(p_test_a_poss.intersection(known_a_numbers)) == 1:
                newly_solved_p.add(p_test)
        known_p_numbers.update(newly_solved_p)
        
    print("They never give an answer.")

solve_puzzle()