import fractions

def solve_cathedral_echo():
    """
    Solves the Cathedral's Echo word problem to find how many pipes the tuner must find.
    """
    # Number of pipes that are still in tune
    in_tune_pipes = 200

    # The fractions of pipes that fell out of tune
    frac_out_1 = fractions.Fraction(1, 3)
    frac_out_2 = fractions.Fraction(2, 5)

    # Step 1: Calculate the total fraction of pipes that are out of tune
    total_frac_out = frac_out_1 + frac_out_2

    # Step 2: Determine the fraction of pipes that are still in tune
    frac_in_tune = 1 - total_frac_out

    # Step 3: Calculate the total number of pipes
    # total_pipes = in_tune_pipes / frac_in_tune
    total_pipes = int(in_tune_pipes / frac_in_tune)

    # Step 4: Find the number of out-of-tune ("lost") pipes
    out_of_tune_pipes = total_pipes - in_tune_pipes

    # Step 5: Calculate the number the tuner must find (half of the lost pipes)
    tuner_must_find = out_of_tune_pipes // 2

    print(f"First, we calculate the total number of pipes in the cathedral.")
    print(f"The total number of pipes is {total_pipes}.")
    print(f"\nNext, we find the number of pipes that are out of tune ('lost').")
    print(f"This is the total number of pipes minus those still in tune: {total_pipes} - {in_tune_pipes} = {out_of_tune_pipes}")
    print(f"\nThe question asks for half of the lost pipes.")
    print(f"The final calculation is:")
    print(f"{out_of_tune_pipes} / 2 = {tuner_must_find}")

solve_cathedral_echo()