import math

def solve_integral():
    """
    Solves the integral by symbolic manipulation and using known mathematical identities.
    """
    
    print("The problem is to evaluate the integral I = integral from 0 to infinity of sum from n=1 to infinity of log(cos(x/2^n)) dx.")
    print("-" * 20)

    # Step 1: Rewrite the sum inside the integral.
    # The sum of logarithms is the logarithm of the product of their arguments.
    print("Step 1: The integrand is a sum of logarithms, which can be expressed as the logarithm of an infinite product:")
    print("  sum_{n=1 to inf} log(cos(x/2^n)) = log( product_{n=1 to inf} cos(x/2^n) )")
    print("-" * 20)

    # Step 2: Use the known identity for the infinite product of cosines.
    print("Step 2: The infinite product has a famous closed-form identity:")
    print("  product_{n=1 to inf} cos(x/2^n) = sin(x) / x")
    print("-" * 20)

    # Step 3: Substitute the result back into the integrand.
    print("Step 3: Substituting the product identity into the expression for the sum gives:")
    print("  Integrand = log(sin(x) / x)")
    print("-" * 20)

    # Step 4: Rewrite the integral.
    # This step assumes that the interchange of summation and integration is permissible.
    print("Step 4: The original integral transforms into:")
    print("  I = integral_{0 to inf} log(sin(x)/x) dx")
    print("-" * 20)

    # Step 5: State the value of the known definite integral.
    # This integral is a standard result known as Lobachevsky's integral.
    print("Step 5: The value of this definite integral is a well-known mathematical result:")
    final_value = -math.pi / 2
    numerator = "-pi"
    denominator = 2
    print(f"  integral_{{0 to inf}} log(sin(x)/x) dx = {numerator} / {denominator}")
    print("-" * 20)

    # Final Answer
    print("Therefore, the value of the integral is:")
    print(f"{numerator} / {denominator} = {final_value}")

if __name__ == '__main__':
    solve_integral()