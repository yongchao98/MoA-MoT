def solve_cathedral_echo():
    """
    Calculates the number of pipes the tuner needs to find based on the poem.
    """
    # Number of pipes that are still in tune.
    pipes_in_tune = 200

    # Fractions of pipes that fell out of tune.
    fraction_out_1 = 1/3
    fraction_out_2 = 2/5

    # Calculate the total fraction of pipes that are out of tune.
    # 1/3 + 2/5 = 5/15 + 6/15 = 11/15
    total_fraction_out = fraction_out_1 + fraction_out_2

    # The remaining fraction of pipes are the ones still in tune.
    # 1 - 11/15 = 4/15
    fraction_in_tune = 1 - total_fraction_out

    # Calculate the total number of pipes. If 4/15 of the total is 200,
    # then the total is 200 / (4/15).
    # 200 * (15/4) = 50 * 15 = 750
    total_pipes = pipes_in_tune / fraction_in_tune

    # Calculate the number of "lost" pipes (total out-of-tune pipes).
    # 750 (total) - 200 (in tune) = 550
    lost_pipes = total_pipes - pipes_in_tune

    # The tuner needs to find half of the lost pipes.
    # 550 / 2 = 275
    pipes_to_find = lost_pipes / 2
    
    # Print the final equation with the calculated numbers as requested.
    # The final step is to find half of the lost pipes.
    print(f"The final calculation is based on the number of 'lost' pipes, which is {int(lost_pipes)}.")
    print(f"The number of pipes the tuner must find is half of the lost pipes:")
    print(f"{int(lost_pipes)} / 2 = {int(pipes_to_find)}")

solve_cathedral_echo()
<<<275>>>