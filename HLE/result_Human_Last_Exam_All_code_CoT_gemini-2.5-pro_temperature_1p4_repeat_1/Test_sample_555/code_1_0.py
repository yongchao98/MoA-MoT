from fractions import Fraction

def solve_dice_problem():
    """
    Calculates the minimal expected number of rolls of a 7-sided die
    to generate a uniform random digit (0-9).
    """
    # Number of sides on the die
    sides = 7
    # Number of target digits
    targets = 10
    
    # We use a strategy based on n=2 rolls.
    n = 2
    total_outcomes = sides ** n

    # After n=2 rolls, we have an outcome v from 0 to 48.
    # The outcome is decisive if floor(10*v/49) == floor(10*(v+1)/49).
    # We count the number of ambiguous outcomes, where this is not true.
    ambiguous_outcomes = 0
    for v in range(total_outcomes):
        if int(targets * v / total_outcomes) != int(targets * (v + 1) / total_outcomes):
            ambiguous_outcomes += 1
            
    successful_outcomes = total_outcomes - ambiguous_outcomes

    prob_success = Fraction(successful_outcomes, total_outcomes)
    prob_failure = Fraction(ambiguous_outcomes, total_outcomes)

    # In case of an ambiguous outcome, we need to solve a subproblem.
    # This subproblem is to generate a Bernoulli(k/10) variable.
    # The expected number of rolls for this, E_B, can be shown to be 7/6
    # for any k from 1 to 9. We can derive this from E_B = 1 + (1/7)*E_B
    # for the p=1/2 case (which has a base-7 expansion of 0.333...).
    E_B = Fraction(7, 6)

    # The total expected number of rolls E is given by the sum of expectations:
    # E = (rolls_on_success * prob_success) + (rolls_on_failure + E_B) * prob_failure
    rolls_spent = n
    E = rolls_spent * prob_success + (rolls_spent + E_B) * prob_failure
    
    # Output the logic and calculation step-by-step.
    print("The optimal strategy involves taking n=2 rolls at a time.")
    print(f"This gives {total_outcomes} total outcomes.")
    print(f"Of these, {successful_outcomes} outcomes are decisive, while {ambiguous_outcomes} are ambiguous.")
    print(f"The probability of a decisive outcome is {successful_outcomes}/{total_outcomes}.")
    print("\nIf the outcome is ambiguous, an average of E_B additional rolls are needed.")
    print(f"This value, E_B, is calculated to be {E_B}.")
    
    print("\nThe total expected number of rolls E is calculated as:")
    print(f"E = (rolls spent) * P(success) + (rolls spent + E_B) * P(failure)")

    # Show the final equation with numbers
    term1_val = rolls_spent * prob_success
    term2_val_factor = rolls_spent + E_B
    
    print("\nSubstituting the values:")
    print(f"E = {rolls_spent} * {prob_success} + ({rolls_spent} + {E_B}) * {prob_failure}")
    print(f"E = {term1_val} + ({term2_val_factor}) * {prob_failure}")
    
    term2_val = term2_val_factor * prob_failure
    print(f"E = {term1_val} + {term2_val}")
    
    # Combine to a single fraction
    final_E = term1_val + term2_val
    print(f"E = {final_E}")
    
    print(f"\nThe minimal expected value of rolls is the simplified fraction: {final_E.numerator}/{final_E.denominator}")

solve_dice_problem()
<<<31/14>>>