import math

def solve_probability():
    """
    Calculates the probability as described in the problem.
    
    The plan is as follows:
    1. Calculate S, the total number of ways to distribute the items.
       S = 25! / (5!)^5
    2. Calculate F, the number of favorable distributions.
       A favorable distribution requires each individual to have a unique specialized item type.
       - There are 5! ways to assign a unique specialized type to each individual.
       - For a fixed assignment (e.g., individual i specializes in type i), we count the number of ways.
         This corresponds to specific configurations (matrices of counts C).
         - Case 1: Each individual i gets 5 of type i. The count matrix is 5*I.
           The number of ways W(C) is 1.
         - Case 2: Each individual i gets 4 of type i, and the remaining items are swapped in a derangement.
           There are D_5 = 44 derangements of 5 items.
           The number of ways W(C) for each is 5^5.
       - Total favorable ways for a fixed assignment is 1 + 44 * 5^5.
       - F = 5! * (1 + 44 * 5^5)
    3. The probability is P = F / S.
    """
    
    # Numerator of S
    s_numerator = math.factorial(25)
    
    # Denominator of S
    s_denominator = math.factorial(5)**5
    
    # Total number of distributions S
    S = s_numerator // s_denominator
    
    # Number of derangements of 5 items
    D5 = 44
    
    # Number of favorable ways for a fixed specialization
    # Case 1: C = 5I, W(C) = 1
    # Case 2: C = 4I + P_derangement, W(C) = 5^5
    fixed_f_ways = 1 + (D5 * (5**5))
    
    # Total number of favorable distributions F
    F = math.factorial(5) * fixed_f_ways
    
    # The probability P = F / S
    # To keep precision, we can use fractions or simplify
    # P = (5! * (1 + 44 * 5^5)) / (25! / 5!^5)
    # P = (5!^6 * (1 + 44 * 5^5)) / 25!
    
    # We will print the components F and S for the final equation.
    
    print(f"The total number of ways to distribute the items is S.")
    print(f"S = 25! / (5!^5) = {S}")
    print()
    print(f"The number of favorable distributions is F.")
    print(f"F = 5! * (1 + 44 * 5^5) = {F}")
    print()
    print(f"The probability is P = F / S.")
    print(f"P = {F} / {S}")

solve_probability()