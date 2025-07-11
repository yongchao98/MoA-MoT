def solve_cathedral_echo():
    """
    Solves the riddle of the Cathedral's Echo to find the number of pipes for the tuner.
    """
    # The core of the problem is setting up an equation for the total number of pipes (T).
    # T = (pipes in tune) + (pipes out of tune)
    # We are given:
    # pipes in tune = 200
    # pipes out of tune = (1/3 * T) + (2/5 * T)
    #
    # So the equation is: T = 200 + (1/3 * T) + (2/5 * T)
    #
    # To solve for T, we rearrange the equation:
    # T - (1/3 * T) - (2/5 * T) = 200
    # T * (1 - 1/3 - 2/5) = 200
    # T * ( (15 - 5 - 6) / 15 ) = 200
    # T * (4 / 15) = 200
    # T = 200 * (15 / 4)

    in_tune_pipes = 200
    
    # Calculate the total number of pipes (T)
    total_pipes = in_tune_pipes * (15 / 4)
    total_pipes = int(total_pipes)
    
    # Calculate the number of "lost" (out-of-tune) pipes
    lost_pipes = total_pipes - in_tune_pipes
    
    # The tuner needs to find half of the lost pipes
    tuner_finds = lost_pipes / 2
    tuner_finds = int(tuner_finds)
    
    # Print the explanation and the final equation with the numbers
    print("Let T be the total number of pipes.")
    print(f"The equation for the total number of pipes is: T = {in_tune_pipes} + (1/3 * T) + (2/5 * T)")
    print(f"Solving this, we find the total number of pipes is {total_pipes}.")
    print(f"The number of 'lost' (out-of-tune) pipes is the total minus the pipes in tune: {total_pipes} - {in_tune_pipes} = {lost_pipes}.")
    print("\nThe question is: How many must the tuner find when just half the lost realign?")
    print("This is half the number of lost pipes.")
    print("\nFinal calculation:")
    print(f"({total_pipes} - {in_tune_pipes}) / 2 = {tuner_finds}")

solve_cathedral_echo()
<<<275>>>