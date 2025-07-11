import math

def solve_optimal_gammas():
    """
    This function calculates and prints the formulas for the optimal two-step
    gradient descent learning rates (gamma_1, gamma_2) for an M-smooth and
    1-strongly convex function.
    """

    print("The optimal choice for the pair of learning rates (gamma_1, gamma_2) is derived from a")
    print("Chebyshev polynomial approximation problem. The formulas below are expressed in terms of M,")
    print("which is the smoothness constant (and condition number, since mu=1).")
    print("-" * 70)

    # The formulas are derived from the roots of the optimal polynomial.
    # The general formulas for m!=1 are:
    # gamma_{a,b} = (4*(M+m) -/+ 2*sqrt(2)*(M-m)) / (M**2 + 6*M*m + m**2)
    # With m=1, M=kappa as per the problem, we get the following:
    
    # Numerator of gamma_1: 4*(M+1) - 2*sqrt(2)*(M-1) = (4 - 2*sqrt(2))*M + (4 + 2*sqrt(2))
    # Numerator of gamma_2: 4*(M+1) + 2*sqrt(2)*(M-1) = (4 + 2*sqrt(2))*M + (4 - 2*sqrt(2))
    # Denominator for both: M^2 + 6*M + 1
    
    # Let's write the polynomial coefficients for gamma = (c1*M + c2) / (M^2 + c3*M + c4)
    c3 = 6.0
    c4 = 1.0
    sqrt2 = math.sqrt(2)

    # For gamma_1
    c1_g1 = 4.0 - 2.0 * sqrt2
    c2_g1 = 4.0 + 2.0 * sqrt2

    # For gamma_2
    c1_g2 = 4.0 + 2.0 * sqrt2
    c2_g2 = 4.0 - 2.0 * sqrt2

    print("\nThe best choice for the pair (gamma_1, gamma_2) is the set containing the following two values:")

    # Print the equation for gamma_1
    print("\nValue 1:")
    print("  gamma_1 = ( (4 - 2*sqrt(2)) * M + (4 + 2*sqrt(2)) ) / ( M^2 + 6*M + 1 )")
    print("\n  The numbers in this equation are:")
    print(f"    Numerator coefficient of M: 4 - 2*sqrt(2) \t= {c1_g1:.15f}")
    print(f"    Numerator constant term:    4 + 2*sqrt(2) \t= {c2_g1:.15f}")
    print(f"    Denominator coefficient of M: \t\t\t= {c3}")
    print(f"    Denominator constant term: \t\t\t\t= {c4}")

    # Print the equation for gamma_2
    print("\nValue 2:")
    print("  gamma_2 = ( (4 + 2*sqrt(2)) * M + (4 - 2*sqrt(2)) ) / ( M^2 + 6*M + 1 )")
    print("\n  The numbers in this equation are:")
    print(f"    Numerator coefficient of M: 4 + 2*sqrt(2) \t= {c1_g2:.15f}")
    print(f"    Numerator constant term:    4 - 2*sqrt(2) \t= {c2_g2:.15f}")
    print(f"    Denominator coefficient of M: \t\t\t= {c3}")
    print(f"    Denominator constant term: \t\t\t\t= {c4}")

    print("\nNote: The order of gamma_1 and gamma_2 does not affect the final error after two steps,")
    print("so the pair simply consists of these two values.")

if __name__ == '__main__':
    solve_optimal_gammas()
