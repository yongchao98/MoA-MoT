def solve_cathedral_riddle():
    """
    This function solves the riddle from "The Cathedral's Echo".
    It calculates the number of pipes the tuner needs to find based on the provided fractions.
    """
    # The number of pipes that are still in tune.
    pipes_in_tune = 200

    # The poem states "One-third of pipes fell out of tune" and "two-fifths caught the rising moon".
    # We interpret this as two separate fractions of the total pipes becoming out of tune.
    fraction_out_of_tune = 1/3 + 2/5

    # The remaining fraction of pipes are the ones that are in tune.
    fraction_in_tune = 1 - fraction_out_of_tune

    # We can now find the total number of pipes.
    # fraction_in_tune * total_pipes = pipes_in_tune
    # total_pipes = pipes_in_tune / fraction_in_tune
    total_pipes = pipes_in_tune / fraction_in_tune

    # The "lost" pipes are all the pipes that are out of tune.
    lost_pipes = total_pipes - pipes_in_tune

    # The tuner must find "half the lost".
    pipes_to_find = lost_pipes / 2

    # The final output needs to show the equation with the numbers.
    # We convert the numbers to integers as they represent whole pipes.
    lost_pipes_int = int(lost_pipes)
    pipes_to_find_int = int(pipes_to_find)

    print(f"First, we determine the number of lost pipes: {lost_pipes_int}.")
    print(f"The tuner must find half of these.")
    print(f"The final calculation is: {lost_pipes_int} / 2 = {pipes_to_find_int}")

solve_cathedral_riddle()
<<<275>>>