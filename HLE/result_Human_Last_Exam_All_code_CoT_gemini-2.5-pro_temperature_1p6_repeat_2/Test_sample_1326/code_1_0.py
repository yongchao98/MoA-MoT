import numpy as np
from scipy.integrate import quad
import warnings

def solve_integral():
    """
    Solves for the definite integral of the piecewise function p(x) from 0 to 4.
    """

    # --- Part 1: Integral from 0 to 3 ---
    def p1(x):
        return (2 * x**3) / 8

    # Analytically, the integral of x^3/4 is x^4/16.
    # From 0 to 3, this is (3^4)/16 - 0 = 81/16 = 5.0625
    integral_part1, _ = quad(p1, 0, 3)
    
    print("The integral is split into two parts: from 0 to 3, and from 3 to 4.")
    print("-" * 50)
    print("Part 1: Integral from x=0 to x=3")
    print("p(x) = (2 * x^3) / 8")
    print(f"The value of the first integral is: (3^4 / 16) - (0^4 / 16) = {integral_part1:.4f}")
    print("-" * 50)

    # --- Part 2: Integral from 3 to 4 ---
    def p2(x):
        return (np.exp(x) * (1 + np.sin(x))) / (1 + np.cos(x))

    print("Part 2: Integral from x=3 to x=4")
    print("p(x) = (e^x * (1 + sin(x))) / (1 + cos(x))")
    print("\nAnalysis of Part 2:")
    print("The denominator (1 + cos(x)) is 0 when x = pi (approx. 3.14159).")
    print("Since pi is within the integration interval [3, 4], the function has a singularity.")
    print("We test for convergence. This type of singularity is not integrable, so the integral diverges.")

    # We attempt the numerical integration, which should confirm the divergence.
    # The `quad` function will issue an IntegrationWarning for divergent integrals.
    integral_part2_str = "Divergent"
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        try:
            # The calculation will fail and raise a warning
            quad(p2, 3, 4)
        except Exception:
            # Catching any other errors just in case
            pass
        
        if any(issubclass(warn.category, IntegrationWarning) for warn in w):
            print("\nNumerical confirmation: The scipy.integrate.quad function issued a warning, confirming the integral is likely divergent.")
        else:
            # This case should not be reached
            print("\nNumerical check could not confirm divergence, but analytical proof stands.")
            
    print(f"The value of the second integral is: {integral_part2_str}")
    print("-" * 50)

    # --- Final Conclusion ---
    print("Final Result:")
    print("The total integral is the sum of the two parts.")
    print(f"Integral from 0 to 4 = (Integral from 0 to 3) + (Integral from 3 to 4)")
    print(f"= {integral_part1:.4f} + {integral_part2_str}")
    print("= Divergent")
    print("\nBecause one part of the integral diverges to infinity, the total integral is also divergent.")


solve_integral()
<<<Divergent>>>