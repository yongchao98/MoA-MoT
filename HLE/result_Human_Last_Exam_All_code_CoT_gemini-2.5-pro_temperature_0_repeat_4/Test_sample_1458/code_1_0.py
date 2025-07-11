def solve_cathedral_echo():
    """
    Solves the word problem from "The Cathedral's Echo".
    """
    # Known values from the poem
    in_tune_pipes = 200
    
    # The fractions of pipes that fell out of tune
    fraction_thunder = 1/3
    fraction_moon = 2/5
    
    # Step 1: Calculate the total fraction of pipes that are out of tune.
    # We add the two fractions: 1/3 + 2/5 = 5/15 + 6/15 = 11/15
    fraction_lost = fraction_thunder + fraction_moon
    
    # Step 2: Calculate the fraction of pipes that are still in tune.
    # This is the remainder: 1 - 11/15 = 4/15
    fraction_in_tune = 1 - fraction_lost
    
    # Step 3: Calculate the total number of pipes.
    # If (4/15) of the total pipes is 200, we can find the total.
    # Total = 200 / (4/15)
    total_pipes = in_tune_pipes / fraction_in_tune
    
    # Step 4: Calculate the number of "lost" (out of tune) pipes.
    lost_pipes = total_pipes - in_tune_pipes
    
    # Step 5: The tuner needs to find half of the lost pipes.
    tuner_finds = lost_pipes / 2
    
    # Print the step-by-step calculation
    print("First, we find the total number of pipes in the cathedral.")
    print(f"The fraction of pipes out of tune is 1/3 + 2/5 = 11/15.")
    print(f"The fraction of pipes still in tune is 1 - 11/15 = 4/15.")
    print(f"Since {in_tune_pipes} pipes represent 4/15 of the total, the total number of pipes is:")
    print(f"Total Pipes = {in_tune_pipes} / (4/15) = {int(total_pipes)}")
    print("")
    print("Next, we find the number of 'lost' pipes (those out of tune).")
    print(f"Lost Pipes = Total Pipes - In-Tune Pipes")
    print(f"Lost Pipes = {int(total_pipes)} - {in_tune_pipes} = {int(lost_pipes)}")
    print("")
    print("Finally, the tuner must find half of the lost pipes.")
    print(f"Pipes to Find = {int(lost_pipes)} / 2 = {int(tuner_finds)}")

solve_cathedral_echo()
<<<275>>>