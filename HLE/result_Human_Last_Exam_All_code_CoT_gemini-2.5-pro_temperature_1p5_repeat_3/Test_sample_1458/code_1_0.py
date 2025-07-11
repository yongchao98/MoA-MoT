import fractions

def solve_cathedral_echo():
    """
    Solves the Cathedral's Echo riddle by calculating the number of pipes a tuner must find.
    """
    # Known values from the poem
    pipes_in_tune = 200
    frac1_out = fractions.Fraction(1, 3)
    frac2_out = fractions.Fraction(2, 5)

    # Step 1: Calculate the total fraction of pipes out of tune.
    total_frac_out = frac1_out + frac2_out
    print(f"First, we calculate the total fraction of pipes that went out of tune.")
    print(f"{frac1_out.numerator}/{frac1_out.denominator} + {frac2_out.numerator}/{frac2_out.denominator} = {total_frac_out.numerator}/{total_frac_out.denominator}")

    # Step 2: Calculate the fraction of pipes still in tune.
    frac_in_tune = 1 - total_frac_out
    print(f"\nThe fraction of pipes still in tune is:")
    print(f"1 - {total_frac_out.numerator}/{total_frac_out.denominator} = {frac_in_tune.numerator}/{frac_in_tune.denominator}")

    # Step 3: Determine the total number of pipes.
    total_pipes = int(pipes_in_tune / frac_in_tune)
    print(f"\nWe know {pipes_in_tune} pipes are in tune, which is {frac_in_tune.numerator}/{frac_in_tune.denominator} of the total. So, the total number of pipes is:")
    print(f"{pipes_in_tune} / ({frac_in_tune.numerator}/{frac_in_tune.denominator}) = {total_pipes}")

    # Step 4: Calculate the total number of "lost" pipes.
    lost_pipes = int(total_pipes * total_frac_out)
    print(f"\nThe number of pipes that went out of tune ('the lost') is:")
    print(f"{total_pipes} * ({total_frac_out.numerator}/{total_frac_out.denominator}) = {lost_pipes}")
    
    # Step 5: Calculate the final number of pipes the tuner must find.
    tuner_must_find = int(lost_pipes / 2)
    print(f"\nThe tuner must find the pipes that remain out of tune after half have been realigned. This is half of the lost pipes:")
    print(f"{lost_pipes} / 2 = {tuner_must_find}")

solve_cathedral_echo()
<<<275>>>