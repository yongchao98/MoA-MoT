def solve_schwartz_function_problem():
    """
    This script explains the solution to the following problem:
    Suppose f: R -> R is a Schwartz class function such that the integral of x^k * f(x)
    over R is 0 for all non-negative integers k. Does it follow that f = 0?
    """

    print("--- The Proof ---")
    print("\nThe answer is YES. If all moments of a Schwartz function are zero, the function must be identically zero.")
    print("Here is the step-by-step proof:\n")

    print("Step 1: The Fourier Transform and its Derivatives")
    print("Let f(x) be a Schwartz function. Its Fourier transform is given by:")
    print("  f_hat(ξ) = ∫ f(x) * e^(-2πixξ) dx")
    print("The k-th derivative of the Fourier transform can be found by differentiating under the integral sign:")
    print("  f_hat^(k)(ξ) = ∫ f(x) * (-2πix)^k * e^(-2πixξ) dx\n")

    print("Step 2: Evaluating Derivatives at the Origin")
    print("Let's evaluate the k-th derivative at ξ = 0:")
    print("  f_hat^(k)(0) = ∫ f(x) * (-2πix)^k * e^(0) dx")
    print("This simplifies to an equation relating the derivative of the transform to the moment of the function:")
    print("  f_hat^(k)(0) = (-2πi)^k * ∫ x^k * f(x) dx\n")

    print("Step 3: Applying the Given Condition")
    print("We are given that all moments of f(x) are zero:")
    print("  ∫ x^k * f(x) dx = 0  (for all integers k ≥ 0)")
    print("Substituting this into our equation from Step 2, we get:")
    # Printing the equation as requested
    print("  f_hat^(k)(0) = (-2πi)^k * 0")
    print("This means all derivatives of f_hat(ξ) at the origin are zero:")
    print("  f_hat^(k)(0) = 0  (for all k ≥ 0)\n")

    print("Step 4: Using the Analytic Property of the Fourier Transform")
    print("A crucial property of Schwartz functions is that their Fourier transforms are not only smooth but also analytic.")
    print("This means f_hat(ξ) can be represented by its Taylor series around ξ = 0:")
    print("  f_hat(ξ) = Σ_{k=0 to ∞} [f_hat^(k)(0) / k!] * ξ^k")
    print("Since every term f_hat^(k)(0) is 0, the entire series is 0:")
    # Printing the equation as requested
    print("  f_hat(ξ) = Σ_{k=0 to ∞} [0 / k!] * ξ^k = 0")
    print("Thus, the Fourier transform f_hat(ξ) is identically zero for all ξ.\n")

    print("Step 5: Conclusion via Inverse Fourier Transform")
    print("The Fourier transform is an invertible mapping on the space of Schwartz functions.")
    print("If the transform f_hat(ξ) is the zero function, then the original function f(x) must also be the zero function, as shown by the inverse transform:")
    print("  f(x) = ∫ f_hat(ξ) * e^(2πixξ) dξ = ∫ 0 * e^(2πixξ) dξ = 0")
    print("Therefore, f(x) = 0 for all x.\n")

if __name__ == '__main__':
    solve_schwartz_function_problem()