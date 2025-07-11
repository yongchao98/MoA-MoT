import collections

def solve_number_puzzle():
    """
    Simulates the logic puzzle step-by-step to find when an answer is given.
    """
    aric_n = 13
    pi_n = 10
    possible_sums = {20, 23}

    # Aric's initial possibilities for Pi's number
    aric_p_poss = {s - aric_n for s in possible_sums if s - aric_n > 0}
    # Pi's initial possibilities for Aric's number
    pi_a_poss = {s - pi_n for s in possible_sums if s - pi_n > 0}

    # These sets will store the numbers that are publicly known to be impossible
    # for Aric or Pi to have, based on their previous passes.
    eliminated_a = set()
    eliminated_p = set()

    print("Initial State:")
    print(f"Aric's secret number is {aric_n}.")
    print(f"Pi's secret number is {pi_n}.")
    print(f"They both know the sum is one of {sorted(list(possible_sums))}.")
    print(f"Aric considers Pi's number could be {sorted(list(aric_p_poss))}.")
    print(f"Pi considers Aric's number could be {sorted(list(pi_a_poss))}.")
    print("-" * 30)

    # We will simulate for a maximum of 10 days
    for day in range(1, 11):
        print(f"--- Day {day} ---")

        # === Aric's Turn ===
        # Aric checks if the accumulated knowledge about Pi's number eliminates one of his possibilities.
        if len(aric_p_poss.intersection(eliminated_p)) == 1:
            print("Aric's Turn: He analyzes the new information.")
            eliminated_possibility = aric_p_poss.intersection(eliminated_p).pop()
            print(f"He knows from previous rounds that Pi cannot have the number {eliminated_possibility}.")
            print(f"This rules out one of his two scenarios (sum = {aric_n + eliminated_possibility}).")
            
            # The surviving possibility gives the answer
            surviving_p = (aric_p_poss - eliminated_p).pop()
            final_sum = aric_n + surviving_p
            print(f"Aric deduces Pi's number must be {surviving_p}.")
            print(f"Therefore, Aric confidently states the sum is {aric_n} + {surviving_p} = {final_sum}.")
            print(f"\nThey give an answer on day {day}.")
            # The puzzle asks for the final number to be wrapped.
            print(f"\n<<<{day}>>>")
            return

        # If Aric cannot solve it, he passes. This pass is new information.
        print(f"Aric's Turn: His possible numbers for Pi, {sorted(list(aric_p_poss))}, are both still plausible.")
        print("Aric passes.")

        # The world learns that Aric's number is not one that would have let him solve it on this turn.
        # We find which numbers those are and add them to the eliminated set.
        newly_eliminated_a = set()
        # A number a' is eliminated if having it would have led to a solution for Aric this turn.
        for a_prime in range(1, max(possible_sums)):
            # Skip numbers already known to be impossible
            if a_prime in eliminated_a:
                continue
            
            # What would this hypothetical Aric's possibilities for Pi's number be?
            p_opts = {s - a_prime for s in possible_sums if s - a_prime > 0}
            if len(p_opts) < 2: continue # He would have solved it on day 1 anyway.
            
            # If exactly one of these possibilities was already eliminated, he would know the answer.
            if len(p_opts.intersection(eliminated_p)) == 1:
                newly_eliminated_a.add(a_prime)
        
        # Add the newly deduced impossible numbers to the public knowledge
        if newly_eliminated_a:
            print(f"Pi learns from Aric's pass that Aric's number cannot be in {sorted(list(newly_eliminated_a))}.")
        eliminated_a.update(newly_eliminated_a)
        print()

        # === Pi's Turn ===
        # Pi checks if the accumulated knowledge about Aric's number eliminates one of her possibilities.
        if len(pi_a_poss.intersection(eliminated_a)) == 1:
            print("Pi's Turn: She analyzes the new information.")
            eliminated_possibility = pi_a_poss.intersection(eliminated_a).pop()
            print(f"She knows from previous rounds that Aric cannot have the number {eliminated_possibility}.")
            print(f"This rules out one of her two scenarios (sum = {pi_n + eliminated_possibility}).")
            
            # The surviving possibility gives the answer
            surviving_a = (pi_a_poss - eliminated_a).pop()
            final_sum = pi_n + surviving_a
            print(f"Pi deduces Aric's number must be {surviving_a}.")
            print(f"Therefore, Pi confidently states the sum is {surviving_a} + {pi_n} = {final_sum}.")
            print(f"\nThey give an answer on day {day}.")
            # The puzzle asks for the final number to be wrapped.
            print(f"\n<<<{day}>>>")
            return

        # If Pi cannot solve it, she passes.
        print(f"Pi's Turn: Her possible numbers for Aric, {sorted(list(pi_a_poss))}, are both still plausible.")
        print("Pi passes.")
        
        # The world learns that Pi's number is not one that would have let her solve it.
        newly_eliminated_p = set()
        for p_prime in range(1, max(possible_sums)):
            if p_prime in eliminated_p:
                continue
            
            a_opts = {s - p_prime for s in possible_sums if s - p_prime > 0}
            if len(a_opts) < 2: continue
                
            if len(a_opts.intersection(eliminated_a)) == 1:
                newly_eliminated_p.add(p_prime)
                
        if newly_eliminated_p:
            print(f"Aric learns from Pi's pass that Pi's number cannot be in {sorted(list(newly_eliminated_p))}.")
        eliminated_p.update(newly_eliminated_p)
        print("-" * 30)

    print("The simulation ended without a resolution.")
    print("\n<<<NEVER>>>")


if __name__ == '__main__':
    solve_number_puzzle()