import fractions

def solve_cathedral_riddle():
    """
    This function solves the word problem from the poem "The Cathedral's Echo".
    It calculates the number of pipes a tuner must find based on the provided numbers and fractions.
    """
    
    # Numbers and fractions from the poem
    in_tune_pipes = 200
    fraction_out_1 = fractions.Fraction(1, 3)
    fraction_out_2 = fractions.Fraction(2, 5)

    # --- Step 1 & 2: Calculate the fraction of pipes in and out of tune ---
    total_fraction_out = fraction_out_1 + fraction_out_2
    fraction_in_tune = 1 - total_fraction_out

    # --- Step 3: Calculate the total number of pipes ---
    # If `fraction_in_tune` of the total pipes is `in_tune_pipes`, then:
    # Total = in_tune_pipes / fraction_in_tune
    total_pipes = in_tune_pipes / fraction_in_tune
    
    # --- Step 4: Calculate the number of out-of-tune pipes ("the lost") ---
    out_of_tune_pipes = total_pipes - in_tune_pipes

    # --- Step 5: Calculate how many pipes the tuner must find (half the lost) ---
    pipes_to_realign = out_of_tune_pipes / 2

    # --- Final Output ---
    # Print the logic and the final equation as requested.
    print("Here is the step-by-step solution to the riddle:")
    
    print(f"\nFirst, we find the total number of out-of-tune pipes.")
    print(f"The total number of out-of-tune pipes ('the lost') can be found by first calculating the total number of pipes.")
    print(f"The fraction of pipes still in tune is 1 - (1/3 + 2/5) = {fraction_in_tune.numerator}/{fraction_in_tune.denominator}.")
    print(f"We know these {int(in_tune_pipes)} pipes represent {fraction_in_tune.numerator}/{fraction_in_tune.denominator} of the total.")
    print(f"Total pipes = {int(in_tune_pipes)} / ({fraction_in_tune.numerator}/{fraction_in_tune.denominator}) = {int(total_pipes)}")
    print(f"The number of lost pipes is {int(total_pipes)} (total) - {int(in_tune_pipes)} (in tune) = {int(out_of_tune_pipes)}.")
    
    print(f"\nThe question is: 'How many must the tuner find when just half the lost realign?'")
    print("This requires calculating half of the out-of-tune pipes.")
    
    print("\nThe final equation is:")
    print(f"{int(out_of_tune_pipes)} / 2 = {int(pipes_to_realign)}")

solve_cathedral_riddle()