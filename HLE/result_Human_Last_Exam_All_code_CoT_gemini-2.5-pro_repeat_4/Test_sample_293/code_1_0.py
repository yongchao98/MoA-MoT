import math

def solve_integral():
    """
    Solves the integral by simplifying the sum and using a known mathematical result.
    """
    # The problem is to evaluate the integral:
    # I = integral from 0 to infinity of sum_{n=1 to infinity} log(cos(x/2^n)) dx

    # Step 1: Simplify the sum inside the integral.
    # The sum is S(x) = sum_{n=1 to inf} log(cos(x/2^n)).
    # Using the logarithm property log(a) + log(b) = log(ab), the sum becomes:
    # S(x) = log( product_{n=1 to inf} cos(x/2^n) )

    # Step 2: Evaluate the infinite product.
    # The infinite product is a known identity:
    # product_{n=1 to inf} cos(x/2^n) = sin(x)/x.
    # This identity is derived from the sine double-angle formula, sin(a) = 2*sin(a/2)*cos(a/2).
    # So, the integrand simplifies to log(sin(x)/x).

    # Step 3: The integral becomes:
    # I = integral from 0 to infinity of log(sin(x)/x) dx.

    # Step 4: Evaluate the definite integral.
    # This is a known, non-trivial integral in mathematical analysis. Its value is -pi/2.
    # The evaluation requires advanced techniques because the integral is conditionally convergent.

    # The final value is -pi / 2.
    pi_value = math.pi
    sign = -1
    numerator = pi_value
    denominator = 2
    
    result = sign * numerator / denominator

    print("The integral is evaluated by first simplifying the sum into log(sin(x)/x).")
    print("The problem then reduces to calculating the integral of log(sin(x)/x) from 0 to infinity.")
    print("The value of this integral is a known result in mathematics.")
    
    print("\nThe final equation is:")
    print(f"Value = ({sign}) * ({numerator}) / ({denominator})")
    print(f"Value = {result}")

solve_integral()