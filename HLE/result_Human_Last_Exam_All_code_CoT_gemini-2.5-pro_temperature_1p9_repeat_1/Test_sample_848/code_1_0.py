import math

def solve():
    """
    This function calculates the integer part of 10^4 * lim_{N->inf} F(N)/ln(N).
    The limit is given by the sum of contributions from two infinite families of solutions.
    """

    # Numbers from the first family of solutions (k=13)
    c1 = 13
    d1 = 165
    
    # Numbers from the second family of solutions (k=5)
    c2 = 5
    d2 = 21

    # Denominators inside the logs
    den = 2
    
    # The constants defining the asymptotic growth rate of the solution sequences
    lambda1 = (c1 + math.sqrt(d1)) / den
    lambda2 = (c2 + math.sqrt(d2)) / den
    
    # The formula for the limit is 2 * (1/ln(lambda1) + 1/ln(lambda2))
    factor1 = 10000
    factor2 = 2

    # The expression derived
    lim_val = factor2 * (1 / math.log(lambda1) + 1 / math.log(lambda2))
    
    # The final value to be computed
    result = factor1 * lim_val

    print("The final expression to calculate is:")
    print(f"floor({factor1} * {factor2} * (1 / ln(({c1} + sqrt({d1})) / {den}) + 1 / ln(({c2} + sqrt({d2})) / {den})))")
    
    print("\nCalculated value:")
    print(math.floor(result))

solve()
