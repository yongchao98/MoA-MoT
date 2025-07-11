import fractions

def solve_cathedral_puzzle():
    """
    This function solves the word puzzle about the cathedral pipes.
    """
    # Number of pipes that are still in tune
    in_tune_pipes = 200

    # Fractions of pipes that went out of tune
    out_of_tune_fraction1 = fractions.Fraction(1, 3)
    out_of_tune_fraction2 = fractions.Fraction(2, 5)

    # Step 1: Calculate the total fraction of pipes that are out of tune
    total_out_of_tune_fraction = out_of_tune_fraction1 + out_of_tune_fraction2
    print(f"First, we find the total fraction of pipes that are out of tune:")
    print(f"1/3 + 2/5 = {total_out_of_tune_fraction.numerator}/{total_out_of_tune_fraction.denominator}\n")

    # Step 2: Calculate the fraction of pipes that are still in tune
    in_tune_fraction = 1 - total_out_of_tune_fraction
    print(f"Next, we find the fraction of pipes that are still in tune:")
    print(f"1 - {total_out_of_tune_fraction.numerator}/{total_out_of_tune_fraction.denominator} = {in_tune_fraction.numerator}/{in_tune_fraction.denominator}\n")

    # Step 3: Calculate the total number of pipes
    # in_tune_fraction * total_pipes = in_tune_pipes
    total_pipes = in_tune_pipes / in_tune_fraction
    print(f"We know {in_tune_pipes} pipes are in tune, which is {in_tune_fraction.numerator}/{in_tune_fraction.denominator} of the total.")
    print(f"So, the total number of pipes is {in_tune_pipes} / ({in_tune_fraction.numerator}/{in_tune_fraction.denominator}) = {int(total_pipes)}.\n")

    # Step 4: Calculate the total number of "lost" (out-of-tune) pipes
    lost_pipes = total_pipes - in_tune_pipes
    print(f"The number of 'lost' pipes is the total minus the in-tune pipes:")
    print(f"{int(total_pipes)} - {in_tune_pipes} = {int(lost_pipes)}\n")

    # Step 5: Answer the final question
    # The tuner must find half of the lost pipes.
    pipes_to_find = lost_pipes / 2
    print(f"The tuner must find half of the lost pipes. The final calculation is:")
    print(f"{int(lost_pipes)} / 2 = {int(pipes_to_find)}")

solve_cathedral_puzzle()
<<<275>>>