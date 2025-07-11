import math
import argparse

def calculate_l(n: int):
    """
    Calculates the value of l(n) based on the derived formula.

    The formula for l(n) is:
    l(n) = (2/n^2) * (n^2 + 1 - (2n-1)*sqrt(n^2-n+1))
    
    This function computes the value and prints the components of the equation.
    """
    if n < 5:
        print("Error: The function l(n) is defined for n >= 5.")
        return

    # Calculate the components of the formula
    n_squared = n * n
    sqrt_term_inner = n_squared - n + 1
    sqrt_term = math.sqrt(sqrt_term_inner)
    
    term1 = n_squared + 1
    term2_factor = 2 * n - 1
    
    numerator = 2 * (term1 - term2_factor * sqrt_term)
    denominator = n_squared
    
    l_n = numerator / denominator

    # Print the results in a structured way
    print(f"Calculating l(n) for n = {n}")
    print("-" * 30)
    print("The formula is: l(n) = (2/n^2) * ( (n^2 + 1) - (2n - 1) * sqrt(n^2 - n + 1) )")
    print("\nIntermediate values:")
    print(f"  n^2             = {n_squared}")
    print(f"  n^2 + 1         = {term1}")
    print(f"  2n - 1          = {term2_factor}")
    print(f"  n^2 - n + 1     = {sqrt_term_inner}")
    print(f"  sqrt(n^2-n+1)   = {sqrt_term}")
    
    print("\nFinal equation with numbers:")
    print(f"  l({n}) = (2 / {denominator}) * ( {term1} - {term2_factor} * {sqrt_term} )")
    
    print("\nFinal result:")
    print(f"  l({n}) = {l_n}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate the exact value of the function l(n).")
    parser.add_argument('n', type=int, help="An integer n where n >= 5.")
    args = parser.parse_args()
    
    calculate_l(args.n)
