def solve_cathedral_riddle():
    """
    Solves the riddle by calculating the number of pipes the tuner must find.
    """
    # The number of pipes that still sing pure (in tune).
    in_tune_pipes = 200

    # The riddle states "One-third of pipes fell out of tune".
    # This means the fraction of pipes still in tune is 1 - 1/3.
    # In fractional terms, this is 2/3.
    fraction_in_tune_numerator = 2
    fraction_in_tune_denominator = 3

    print("The plan is to first find the total number of pipes, then the number of 'lost' pipes, and finally the number the tuner must find.")
    print("-" * 20)
    
    # We can set up the equation: (2/3) * Total_Pipes = 200
    print(f"Step 1: Find the total number of pipes (T).")
    print(f"The {in_tune_pipes} pure pipes represent {fraction_in_tune_numerator}/{fraction_in_tune_denominator} of the total.")
    print(f"The equation is: ({fraction_in_tune_numerator}/{fraction_in_tune_denominator}) * T = {in_tune_pipes}")

    # To find the total number of pipes, we solve for T.
    # T = 200 * (3/2)
    total_pipes = in_tune_pipes * fraction_in_tune_denominator / fraction_in_tune_numerator
    print(f"Solving for T: T = {in_tune_pipes} * {fraction_in_tune_denominator} / {fraction_in_tune_numerator}")
    print(f"Total pipes = {int(total_pipes)}")
    print("-" * 20)

    # The number of "lost" pipes is the total minus those that are in tune.
    lost_pipes = total_pipes - in_tune_pipes
    print(f"Step 2: Find the number of 'lost' (out of tune) pipes.")
    print(f"The number of lost pipes = Total Pipes - In-Tune Pipes")
    print(f"Calculation: {int(total_pipes)} - {in_tune_pipes} = {int(lost_pipes)}")
    print("-" * 20)

    # The question asks how many the tuner finds when "just half the lost realign".
    tuner_finds = lost_pipes / 2
    print("Step 3: Answer the final question.")
    print("How many must the tuner find when just half the lost realign?")
    print(f"The final calculation is: {int(lost_pipes)} / 2 = {int(tuner_finds)}")

solve_cathedral_riddle()