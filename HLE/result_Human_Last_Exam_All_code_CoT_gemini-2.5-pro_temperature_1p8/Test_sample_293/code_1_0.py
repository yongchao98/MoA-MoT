import math

def solve_integral():
    """
    This function explains the step-by-step solution to the integral
    I = integral from 0 to infinity of sum from n=1 to infinity of log(cos(x/2^n)) dx
    and prints the final answer.
    """

    print("Step 1: The integral to be evaluated is:")
    print("I = integral_{0}^{infinity} sum_{n=1}^{infinity} log(cos(x / 2^n)) dx")
    print("-" * 30)

    print("Step 2: We simplify the integrand by converting the sum of logarithms into a logarithm of a product.")
    print("sum_{n=1}^{infinity} log(cos(x / 2^n)) = log(product_{n=1}^{infinity} cos(x / 2^n))")
    print("-" * 30)

    print("Step 3: We evaluate the infinite product using the identity based on sin(2a) = 2sin(a)cos(a).")
    print("The product telescopes and evaluates to sin(x)/x.")
    print("product_{n=1}^{infinity} cos(x / 2^n) = sin(x) / x")
    print("-" * 30)
    
    print("Step 4: Substituting this back, the integral becomes:")
    print("I = integral_{0}^{infinity} log(sin(x) / x) dx")
    print("This is a known result in mathematics (Malmsten's or Lobachevsky's integral).")
    print("-" * 30)
    
    print("Step 5: The value of this integral is -pi / 2.")
    value = -math.pi / 2
    
    print("\nFinal Equation:")
    print(f"integral_{0}^{infinity} log(sin(x) / x) dx = -pi / 2")
    
    print("\nNumerical Value:")
    print(f"The value of the integral is {value}")

if __name__ == "__main__":
    solve_integral()