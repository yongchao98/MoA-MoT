import math

def solve_simplices_count():
    """
    Calculates the number of n-simplices for the given mathematical problem.
    """
    N = 200
    k = 13

    # The number of objects in the category from which we are choosing is m.
    # These are the integers in the range [-N, k].
    m = k - (-N) + 1

    print(f"Given N = {N} and k = {k}, the number of integers to choose from is {k} - (-{N}) + 1 = {m}.")
    print("The number of n-simplices is given by the formula C(m + n, n + 1).\n")

    # We need to calculate for n from 0 to 5.
    for n in range(6):
        # The number of n-simplices is the number of ways to choose n+1 items
        # from m with replacement. This is given by the binomial coefficient:
        # C(m + (n+1) - 1, n+1) which simplifies to C(m + n, n+1).
        
        upper = m + n
        lower = n + 1
        
        # Use math.comb for accurate calculation of combinations.
        result = math.comb(upper, lower)
        
        print(f"For n = {n}:")
        
        # Build the string for the full equation as requested.
        # C(upper, lower) = upper * (upper-1) * ... * (upper-lower+1) / lower!
        numerator_terms = [str(i) for i in range(upper, upper - lower, -1)]
        denominator_terms = [str(i) for i in range(lower, 0, -1)]
        
        numerator_str = " * ".join(numerator_terms)
        denominator_str = " * ".join(denominator_terms)
        
        # The factorial part in the denominator is lower!
        denominator_val = math.factorial(lower)

        if lower > 1:
            print(f"Number of {n}-simplices = C({upper}, {lower}) = ({numerator_str}) / ({denominator_str}) = {result}")
        else:
            # Handle the case for n=0 where the denominator is 1.
            print(f"Number of {n}-simplices = C({upper}, {lower}) = {upper} / {lower} = {result}")
        print("-" * 20)

solve_simplices_count()