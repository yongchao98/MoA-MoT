import math

def solve():
    """
    Calculates the number of n-simplices for the given mathematical construction.
    """
    N = 200
    k = 13
    
    print(f"Given N = {N} and k = {k}, the number of n-simplices is calculated for n = 0 to 5.")
    
    for n in range(6):
        # The formula for the number of n-simplices is C(N - k + n + 1, n + 1)
        # With N=200 and k=13, this is C(188 + n, n + 1)
        
        numerator = 188 + n
        denominator = n + 1
        
        # Using math.comb to calculate the binomial coefficient C(numerator, denominator)
        num_simplices = math.comb(numerator, denominator)
        
        print(f"For n = {n}, the number of {n}-simplices is C({numerator}, {denominator}) = {num_simplices}")

solve()