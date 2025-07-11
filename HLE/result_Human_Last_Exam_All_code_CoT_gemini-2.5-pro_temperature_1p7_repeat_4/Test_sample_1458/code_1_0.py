def solve_cathedral_echo():
    """
    Calculates the number of pipes a tuner must find based on the poem.
    """
    # Known values from the poem
    pipes_in_tune = 200
    fraction_out_thunder = 1/3
    fraction_out_moon = 2/5

    # --- Step 1 & 2: Calculate the total number of pipes ---
    
    # The total fraction of pipes that are out of tune
    total_fraction_out = fraction_out_thunder + fraction_out_moon
    
    # The fraction of pipes that remain in tune is 1 minus the out-of-tune fraction
    fraction_in_tune = 1 - total_fraction_out
    
    # Calculate the total number of pipes using the known number of in-tune pipes
    # total_pipes * fraction_in_tune = pipes_in_tune
    total_pipes = pipes_in_tune / fraction_in_tune
    
    # Ensure the result is an integer, as we are counting pipes
    total_pipes = int(total_pipes)

    # --- Step 3: Calculate the number of lost pipes ---
    lost_pipes = total_pipes - pipes_in_tune

    # --- Step 4: Calculate the final answer ---
    # The tuner needs to find half of the lost pipes
    pipes_to_find = lost_pipes / 2
    pipes_to_find = int(pipes_to_find)

    # --- Output the results and the final equation ---
    print(f"First, we find the total number of 'lost' pipes (those out of tune).")
    print(f"The fraction of pipes remaining in tune is 1 - (1/3 + 2/5) = 4/15.")
    print(f"If 200 pipes represent 4/15 of the total, the total number of pipes is {total_pipes}.")
    print(f"The number of lost pipes is the total minus those in tune: {total_pipes} - {pipes_in_tune} = {lost_pipes}.")
    
    print("\nThe question asks how many pipes the tuner must find, which is half the lost ones.")
    print("The final equation is:")
    print(f"{lost_pipes} / 2 = {pipes_to_find}")


solve_cathedral_echo()
<<<275>>>