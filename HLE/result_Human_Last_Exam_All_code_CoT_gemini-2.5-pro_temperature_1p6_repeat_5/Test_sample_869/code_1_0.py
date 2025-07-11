import math

def solve_probability():
    """
    Calculates the total and favorable number of distributions to find the probability.
    """
    
    # Step 1: Calculate S, the total number of ways to distribute the items.
    # S = 25! / (5!)^5
    try:
        s_numerator = math.factorial(25)
        s_denominator = math.factorial(5)**5
        
        # Ensure integer division
        if s_numerator % s_denominator != 0:
            print("Error: S is not an integer, which should not happen in this problem.")
            return
            
        S = s_numerator // s_denominator
    except OverflowError:
        print("Error: Calculation involves numbers too large to handle.")
        return

    # Step 2: Calculate F, the number of favorable distributions.
    # Based on the simplifying assumption that only "pure" distributions are counted,
    # where each individual gets 5 items of a single, unique type.
    # There are 5! ways to assign the 5 types to the 5 individuals.
    F = math.factorial(5)
    
    # Step 3: Print the final equation for the probability P = F / S.
    print("The total number of ways to distribute the items is S.")
    print(f"S = 25! / (5! * 5! * 5! * 5! * 5!)")
    print(f"S = {s_numerator} / {s_denominator}")
    print(f"S = {S}")
    print()
    print("The number of favorable distributions (under a simplifying assumption) is F.")
    print(f"F = 5!")
    print(f"F = {F}")
    print()
    print("The probability P is F / S.")
    print(f"P = {F} / {S}")

solve_probability()