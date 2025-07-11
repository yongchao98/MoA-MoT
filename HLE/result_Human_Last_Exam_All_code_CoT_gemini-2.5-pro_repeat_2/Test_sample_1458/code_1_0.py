def solve_cathedral_echo():
    """
    Solves the Cathedral's Echo riddle by calculating the number of pipes a tuner must find.
    """
    # Known values from the poem
    in_tune_count = 200
    fraction_out_thunder = 1/3
    fraction_out_moon = 2/5

    # Step 1: Calculate the total fraction of pipes out of tune.
    # The problem states 1/3 and 2/5 of the pipes went out of tune.
    total_fraction_out = fraction_out_thunder + fraction_out_moon
    
    # Step 2: Calculate the fraction of pipes still in tune.
    # This is the remainder of the total pipes.
    fraction_in_tune = 1 - total_fraction_out
    
    # Step 3: Calculate the total number of pipes in the cathedral.
    # We know the count of in-tune pipes (200) represents the fraction calculated above.
    total_pipes = in_tune_count / fraction_in_tune
    
    # Step 4: Calculate the total number of "lost" (out-of-tune) pipes.
    lost_pipes = total_pipes - in_tune_count
    
    # Step 5: Calculate the final answer.
    # The question asks for the number to find to realign half of the lost pipes.
    pipes_to_realign = lost_pipes / 2

    # Print the step-by-step solution
    print("Step 1: Calculate the fraction of pipes out of tune.")
    print(f"Fraction = {int(fraction_out_thunder * 15)}/15 (from thunder) + {int(fraction_out_moon * 15)}/15 (from moon) = {int(total_fraction_out * 15)}/15")
    print("-" * 30)

    print("Step 2: Calculate the fraction of pipes still in tune.")
    print(f"Fraction = 1 - {int(total_fraction_out * 15)}/15 = {int(fraction_in_tune * 15)}/15")
    print("-" * 30)

    print("Step 3: Calculate the total number of pipes.")
    print(f"Total Pipes = {in_tune_count} (in tune) / ({int(fraction_in_tune * 15)}/15) = {int(total_pipes)}")
    print("-" * 30)
    
    print("Step 4: Calculate the number of 'lost' pipes.")
    print(f"Lost Pipes = {int(total_pipes)} (total) - {in_tune_count} (in tune) = {int(lost_pipes)}")
    print("-" * 30)

    print("Step 5: Find the number to realign (half the lost).")
    print(f"Final calculation: {int(lost_pipes)} / 2 = {int(pipes_to_realign)}")
    print("-" * 30)

solve_cathedral_echo()
<<<275>>>