import fractions

def solve_cathedral_riddle():
    """
    Solves the Cathedral's Echo riddle by calculating the number of pipes a tuner must find.
    """
    # Fractions of pipes that went out of tune
    fraction_out_of_tune_1 = fractions.Fraction(1, 3)
    fraction_out_of_tune_2 = fractions.Fraction(2, 5)

    # Total fraction of out-of-tune pipes
    total_fraction_out_of_tune = fraction_out_of_tune_1 + fraction_out_of_tune_2
    print(f"First, we calculate the total fraction of pipes that are out of tune.")
    print(f"This is the sum of {fraction_out_of_tune_1.numerator}/{fraction_out_of_tune_1.denominator} and {fraction_out_of_tune_2.numerator}/{fraction_out_of_tune_2.denominator}, which equals {total_fraction_out_of_tune.numerator}/{total_fraction_out_of_tune.denominator}.")
    print("-" * 20)

    # Fraction of pipes that are still in tune
    fraction_in_tune = 1 - total_fraction_out_of_tune
    print(f"The fraction of pipes still in tune is 1 - {total_fraction_out_of_tune.numerator}/{total_fraction_out_of_tune.denominator}, which is {fraction_in_tune.numerator}/{fraction_in_tune.denominator}.")
    print("-" * 20)
    
    # Number of pipes that are still in tune
    in_tune_pipes_count = 200
    print(f"We are told that {in_tune_pipes_count} pipes still sing pure.")
    
    # Calculate the total number of pipes
    # (fraction_in_tune) * Total = in_tune_pipes_count
    # Total = in_tune_pipes_count / fraction_in_tune
    total_pipes_count = int(in_tune_pipes_count / fraction_in_tune)
    print(f"To find the total number of pipes, we solve: ({fraction_in_tune.numerator}/{fraction_in_tune.denominator}) * Total = {in_tune_pipes_count}.")
    print(f"Total pipes = {in_tune_pipes_count} * {fraction_in_tune.denominator} / {fraction_in_tune.numerator} = {total_pipes_count}.")
    print("-" * 20)

    # Calculate the number of "lost" (out-of-tune) pipes
    lost_pipes_count = total_pipes_count - in_tune_pipes_count
    print(f"The number of 'lost' pipes is the total number of pipes minus the in-tune pipes.")
    print(f"So, {total_pipes_count} - {in_tune_pipes_count} = {lost_pipes_count} lost pipes.")
    print("-" * 20)

    # The final question asks for half of the lost pipes
    pipes_to_realign = lost_pipes_count / 2
    print(f"The final question asks how many pipes the tuner must find to realign half the lost ones.")
    print(f"The final answer is {lost_pipes_count} / 2 = {int(pipes_to_realign)}.")

solve_cathedral_riddle()