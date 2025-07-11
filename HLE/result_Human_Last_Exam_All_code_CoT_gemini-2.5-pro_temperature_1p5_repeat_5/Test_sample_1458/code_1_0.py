def solve_cathedral_echo():
    """
    Solves the Cathedral's Echo word problem by calculating the number of pipes a tuner must find.
    """
    # Known values from the problem
    in_tune_pipes = 200
    fraction_fell_out = 1/3
    fraction_caught_moon = 2/5

    # Step 1: Calculate the fraction of pipes that are still in tune.
    # The total fraction of out-of-tune pipes is the sum of the two groups.
    # We find a common denominator (15) for 3 and 5.
    # 1/3 is 5/15, and 2/5 is 6/15.
    fraction_out_of_tune = fraction_fell_out + fraction_caught_moon  # This is 11/15
    fraction_in_tune = 1 - fraction_out_of_tune # This is 4/15

    # Step 2: Calculate the total number of pipes (T).
    # The equation is: fraction_in_tune * T = in_tune_pipes
    # Or: (4/15) * T = 200
    print("First, we determine the total number of pipes in the organ.")
    print(f"The equation for the pipes still in tune is: (1 - 1/3 - 2/5) * Total_Pipes = {in_tune_pipes}")
    print(f"This simplifies to: (4/15) * Total_Pipes = {in_tune_pipes}")
    total_pipes = in_tune_pipes / fraction_in_tune
    print(f"Solving for the total: Total_Pipes = {in_tune_pipes} * 15 / 4 = {int(total_pipes)}")
    print("-" * 30)

    # Step 3: Calculate the number of "lost" (out-of-tune) pipes.
    lost_pipes = total_pipes - in_tune_pipes
    print("Next, we find the number of 'lost' pipes that are out of tune.")
    print(f"Lost Pipes = Total_Pipes - Pipes_In_Tune")
    print(f"The calculation is: {int(total_pipes)} - {in_tune_pipes} = {int(lost_pipes)}")
    print("-" * 30)

    # Step 4: Calculate the final answer.
    # The tuner must find half of the lost pipes.
    pipes_to_realign = lost_pipes / 2
    print("Finally, we calculate how many pipes the tuner must find, which is half the lost pipes.")
    print(f"The final equation is: {int(lost_pipes)} / 2 = {int(pipes_to_realign)}")

solve_cathedral_echo()
<<<275>>>