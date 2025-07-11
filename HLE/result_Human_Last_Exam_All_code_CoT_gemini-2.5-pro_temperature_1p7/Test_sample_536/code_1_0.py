import math
import numpy as np

def solve_for_inner_product():
    """
    This function solves for the limiting inner product <b_p, z_p> based on the problem's given values.
    """
    # Given limiting values from the problem statement
    A = 0.9375  # lim <h_p, b_p>
    B = 0.9     # lim <h_p, z_p>

    print(f"Let A = lim <h_p, b_p> = {A}")
    print(f"Let B = lim <h_p, z_p> = {B}")
    print(f"We want to find x = lim <b_p, z_p>.")

    # Based on the assumption that the limiting vectors h, b, z are coplanar,
    # we can derive the following quadratic equation for x:
    # x^2 - (2*A*B)x + (A^2 + B^2 - 1) = 0

    # Coefficients of the quadratic equation ax^2 + bx + c = 0
    a_coeff = 1
    b_coeff = -2 * A * B
    c_coeff = A**2 + B**2 - 1

    print("\nDerived quadratic equation for x:")
    # We output each number in the final equation as requested.
    print(f"{a_coeff} * x^2 + ({b_coeff}) * x + ({c_coeff}) = 0")
    
    # Use Fractions for precision to show the numbers from the manual derivation
    from fractions import Fraction
    A_f = Fraction(15, 16)
    B_f = Fraction(9, 10)
    b_coeff_f = -2 * A_f * B_f
    c_coeff_f = A_f**2 + B_f**2 - 1
    print("\nUsing fractions, the equation is:")
    print(f"x^2 - (27/16)x + 4409/6400 = 0")
    print("Multiplying by 6400, we get: 6400x^2 - 10800x + 4409 = 0")
    
    # Calculate the discriminant
    discriminant = b_coeff**2 - 4 * a_coeff * c_coeff

    # Solve for the two roots of x
    root1 = (-b_coeff + math.sqrt(discriminant)) / (2 * a_coeff)
    root2 = (-b_coeff - math.sqrt(discriminant)) / (2 * a_coeff)

    print(f"\nThe two possible solutions for x are:")
    print(f"x1 = {root1:.8f}")
    print(f"x2 = {root2:.8f}")

    # Selection criteria: For h to be geometrically between b and z, the coefficients
    # c1 and c2 in h = c1*b + c2*z must be positive.
    # c1 > 0 implies A - Bx > 0 => x < A/B
    # c2 > 0 implies B - Ax > 0 => x < B/A
    
    thresh1 = A / B
    thresh2 = B / A
    
    print("\nSelection Criteria:")
    print(f"For a physically meaningful solution, we need x < A/B ({thresh1:.4f}) and x < B/A ({thresh2:.4f}).")
    print(f"The stricter condition is x < {min(thresh1, thresh2):.4f}.")
    
    final_answer = None
    # Check which root satisfies the condition
    if root1 < min(thresh1, thresh2):
        final_answer = root1
        print(f"\nRoot x1 = {root1:.8f} satisfies the condition.")
    elif root2 < min(thresh1, thresh2):
        final_answer = root2
        print(f"\nRoot x2 = {root2:.8f} satisfies the condition.")

    if final_answer is not None:
        print(f"\nThe final answer is {final_answer:.8f}")
        # Find the fractional representation
        final_answer_f = Fraction(final_answer).limit_denominator(10000)
        print(f"As a fraction, the answer is: {final_answer_f}")
    else:
        print("\nNeither root satisfies the condition, there might be an issue with the model.")

solve_for_inner_product()