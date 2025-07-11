def solve_cathedral_riddle():
    """
    This script solves the riddle of the Cathedral's Echo by calculating
    the number of pipes a tuner must find.
    """

    # The number of pipes that are still in tune is given.
    in_tune_pipes = 200

    # The riddle implies two groups of pipes went out of tune.
    # The 'one-third' group is clear. The 'two-fifths' is a red herring.
    # The 'one-fourth' mentioned later is the second group.
    # We calculate the total fraction of pipes that are out of tune.
    # Fraction from thunder: 1/3
    # Fraction descending minor scales: 1/4
    out_of_tune_fraction_numerator = 1 * 4 + 1 * 3
    out_of_tune_fraction_denominator = 3 * 4
    
    print("Step 1: Determine the fraction of pipes that are out of tune.")
    print("The total fraction of out-of-tune pipes is the sum of two groups: 1/3 and 1/4.")
    print(f"Fraction out of tune = 1/3 + 1/4 = {out_of_tune_fraction_numerator}/{out_of_tune_fraction_denominator}")
    print("-" * 30)

    # The remaining fraction of pipes must be in tune.
    in_tune_fraction_numerator = out_of_tune_fraction_denominator - out_of_tune_fraction_numerator
    in_tune_fraction_denominator = out_of_tune_fraction_denominator

    print("Step 2: Determine the fraction of pipes that are in tune.")
    print(f"Fraction in tune = 1 - {out_of_tune_fraction_numerator}/{out_of_tune_fraction_denominator} = {in_tune_fraction_numerator}/{in_tune_fraction_denominator}")
    print("-" * 30)

    # We can now calculate the total number of pipes.
    total_pipes = int(in_tune_pipes * in_tune_fraction_denominator / in_tune_fraction_numerator)

    print("Step 3: Calculate the total number of pipes in the cathedral.")
    print(f"We know that {in_tune_fraction_numerator}/{in_tune_fraction_denominator} of the total pipes equals {in_tune_pipes}.")
    print(f"Total pipes = {in_tune_pipes} / ({in_tune_fraction_numerator}/{in_tune_fraction_denominator}) = {total_pipes}")
    print("-" * 30)
    
    # Calculate the number of out-of-tune ("lost") pipes.
    lost_pipes = total_pipes - in_tune_pipes

    print("Step 4: Calculate the total number of 'lost' (out-of-tune) pipes.")
    print(f"Lost pipes = Total pipes - In-tune pipes")
    print(f"Lost pipes = {total_pipes} - {in_tune_pipes} = {lost_pipes}")
    print("-" * 30)
    
    # The final question asks how many pipes the tuner must find, which is half of the lost pipes.
    tuner_finds = int(lost_pipes / 2)
    
    print("Step 5: Calculate the final answer.")
    print("The tuner must find half of the lost pipes.")
    print(f"Pipes to find = {lost_pipes} / 2 = {tuner_finds}")

solve_cathedral_riddle()
<<<140>>>