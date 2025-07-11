def solve_cathedral_echo():
    """
    Calculates the number of pipes the tuner must find based on the riddle.
    """
    # Number of pipes that are confirmed to be in tune.
    pipes_in_tune = 200

    # Fractions of pipes that went out of tune.
    frac1_out = 1/3
    frac2_out = 2/5

    # Step 1: Calculate the total fraction of pipes out of tune.
    # The common denominator for 3 and 5 is 15.
    # 1/3 is 5/15. 2/5 is 6/15.
    total_frac_out = frac1_out + frac2_out
    
    # Step 2: Calculate the fraction of pipes remaining in tune.
    frac_in_tune = 1 - total_frac_out

    # Step 3: Calculate the total number of pipes in the cathedral.
    # We know that the number of pipes in tune (200) is equal to
    # the total number of pipes multiplied by the fraction in tune.
    # So, Total Pipes = Pipes in Tune / Fraction in Tune
    total_pipes = pipes_in_tune / frac_in_tune
    
    # Step 4: Calculate the total number of "lost" (out of tune) pipes.
    lost_pipes = total_pipes * total_frac_out

    # Step 5: Calculate how many pipes the tuner must find.
    # This is half of the lost pipes.
    tuner_must_find = lost_pipes / 2

    # --- Output the step-by-step calculation ---
    print("Let's solve the riddle of the Cathedral's Echo.")
    print("\nFirst, we find the total number of pipes:")
    print(f"Fraction of pipes out of tune = 1/3 + 2/5 = {int(frac1_out*15)}/15 + {int(frac2_out*15)}/15 = {int(total_frac_out*15)}/15")
    print(f"Fraction of pipes in tune = 1 - {int(total_frac_out*15)}/15 = {int(frac_in_tune*15)}/15")
    print(f"We know {int(frac_in_tune*15)}/15 of the pipes equals {pipes_in_tune}.")
    print(f"So, the total number of pipes is {pipes_in_tune} / ({int(frac_in_tune*15)}/15) = {int(total_pipes)}")

    print("\nNext, we find the number of lost pipes:")
    print(f"Number of lost pipes = {int(total_pipes)} (total) - {pipes_in_tune} (in tune) = {int(lost_pipes)}")
    
    print("\nFinally, we calculate how many pipes the tuner must find:")
    print("The tuner must realign half of the lost pipes.")
    print(f"The final equation is: {int(lost_pipes)} / 2 = {int(tuner_must_find)}")


solve_cathedral_echo()
<<<275>>>