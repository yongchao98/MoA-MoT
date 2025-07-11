import numpy as np

def solve_for_q0():
    """
    This function solves the problem by finding the root of the derived polynomial
    and then calculating the final required value.
    """
    # The probability of Alice winning, p, is a function of y = q^4.
    # We set p = 0.95 and solve for y. This leads to a polynomial equation.
    # The equation is: 361*y^3 - 1140*y^2 + 1200*y - 400 = 0
    
    # Coefficients of the polynomial P(y) = 361*y^3 - 1140*y^2 + 1200*y - 400
    coeffs = [361, -1140, 1200, -400]
    
    print("The polynomial equation for y = q^4 is:")
    print(f"{coeffs[0]}*y^3 + ({coeffs[1]})*y^2 + {coeffs[2]}*y + ({coeffs[3]}) = 0")
    
    # Find the roots of the polynomial
    roots = np.roots(coeffs)
    
    # Alice has a non-zero chance of winning only if y > 1/3.
    # The probability p is between 0 and 1, which corresponds to y being between 1/3 and 1.
    # We need to find the real root in this range.
    y0 = -1
    for r in roots:
        if np.isreal(r):
            r_real = np.real(r)
            if (1/3) < r_real < 1:
                y0 = r_real
                break
    
    if y0 == -1:
        print("Could not find a valid root for y0.")
        return

    # We have found y0, which is q0^4. Now we find q0.
    q0 = y0**(1/4)
    
    # The problem asks for floor(100 * q0)
    result = np.floor(100 * q0)
    
    print(f"\nThe relevant root of the equation is y0 = {y0:.6f}")
    print(f"This gives q0 = y0^(1/4) = {q0:.6f}")
    print(f"The value of 100 * q0 is {100 * q0:.6f}")
    print(f"\nThe final answer is floor(100 * q0):")
    print(int(result))

solve_for_q0()