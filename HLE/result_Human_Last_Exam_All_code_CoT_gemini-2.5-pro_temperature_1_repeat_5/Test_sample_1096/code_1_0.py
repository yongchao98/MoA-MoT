import math

def solve_step_sizes():
    """
    This function calculates and prints the formulas for the optimal
    two-step gradient descent step sizes (gamma_1, gamma_2).
    """

    # The problem of finding the best (gamma_1, gamma_2) pair reduces to
    # finding a polynomial P(lambda) = 1 - (g1+g2)*lambda + g1*g2*lambda^2
    # that has the minimum maximum absolute value on the interval [1, kappa],
    # under the constraint P(0)=1.
    # The solution is a scaled Chebyshev polynomial. From its coefficients, we can
    # derive the values for the sum (S) and product (P) of the gammas:
    # S = gamma_1 + gamma_2 = 8*(kappa+1) / (kappa^2 + 6*kappa + 1)
    # P = gamma_1 * gamma_2 = 8 / (kappa^2 + 6*kappa + 1)
    #
    # The individual gammas are the roots of the quadratic equation:
    # z^2 - S*z + P = 0
    # The roots are z = (S +/- sqrt(S^2 - 4P)) / 2.
    # After algebraic simplification, we get the formulas below.

    # Coefficients for the numerators in the expression for gamma_1 and gamma_2
    # gamma_1 = (c1_k * kappa + c1_c) / (kappa^2 + 6*kappa + 1)
    # gamma_2 = (c2_k * kappa + c2_c) / (kappa^2 + 6*kappa + 1)

    c1_k = 4 - 2 * math.sqrt(2)
    c1_c = 4 + 2 * math.sqrt(2)
    c2_k = 4 + 2 * math.sqrt(2)
    c2_c = 4 - 2 * math.sqrt(2)

    # Coefficients for the denominator: kappa^2 + 6*kappa + 1
    d_k2 = 1
    d_k = 6
    d_c = 1

    print("The best choice for the pair (gamma_1, gamma_2) is given by the following formulas,")
    print("where kappa is the condition number M/m (with m=1):")
    print("-" * 60)

    # Outputting each number in the final equation
    print("\nFormula for gamma_1:")
    print("gamma_1 = ((4 - 2*sqrt(2))*kappa + (4 + 2*sqrt(2))) / (kappa^2 + 6*kappa + 1)")
    print(f"Numerical coefficients: gamma_1 = ({c1_k:.6f} * kappa + {c1_c:.6f}) / ({d_k2}*kappa^2 + {d_k}*kappa + {d_c})")

    print("\nFormula for gamma_2:")
    print("gamma_2 = ((4 + 2*sqrt(2))*kappa + (4 - 2*sqrt(2))) / (kappa^2 + 6*kappa + 1)")
    print(f"Numerical coefficients: gamma_2 = ({c2_k:.6f} * kappa + {c2_c:.6f}) / ({d_k2}*kappa^2 + {d_k}*kappa + {d_c})")
    print("-" * 60)
    print("\nNote: The order of gamma_1 and gamma_2 can be swapped as the final result for x_2 is symmetric with respect to them.")

solve_step_sizes()