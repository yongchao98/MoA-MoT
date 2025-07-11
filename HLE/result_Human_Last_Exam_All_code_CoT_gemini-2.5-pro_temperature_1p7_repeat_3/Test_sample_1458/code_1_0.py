import fractions

def solve_cathedral_riddle():
    """
    Solves the word problem from the poem "The Cathedral's Echo".
    """
    # Known values from the poem
    pipes_in_tune = 200
    fraction_out_1 = fractions.Fraction(1, 3)
    fraction_out_2 = fractions.Fraction(2, 5)

    # Step 1: Calculate the total fraction of pipes that are out of tune.
    total_fraction_out_of_tune = fraction_out_1 + fraction_out_2

    # Step 2: Calculate the fraction of pipes that are still in tune.
    fraction_in_tune = 1 - total_fraction_out_of_tune

    # Step 3: Calculate the total number of pipes.
    # (fraction_in_tune * total_pipes) = pipes_in_tune
    # total_pipes = pipes_in_tune / fraction_in_tune
    total_pipes = pipes_in_tune / fraction_in_tune

    # Step 4: Calculate the number of "lost" (out-of-tune) pipes.
    lost_pipes = total_pipes - pipes_in_tune

    # Step 5: The tuner must find half of the lost pipes.
    tuner_finds = lost_pipes / 2

    # Print the final equation as requested.
    # The equation is based on the final step of the problem.
    print("The final calculation to find how many pipes the tuner must realign:")
    print(f"{int(lost_pipes)} / 2 = {int(tuner_finds)}")

solve_cathedral_riddle()