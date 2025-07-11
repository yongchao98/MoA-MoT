import numpy as np
from sympy import S

def solve_problem():
    """
    This function calculates the value of (1 - max|Φ|) based on the problem description.
    """

    # Step 1 & 2: Calculate v1^2 and v2^2 for KNbO3
    c_33 = S(3) * 10**11
    d_33 = S(3) * 10**-11
    k_33 = S(300)
    
    c_44 = S('3/4') * 10**11
    d_24 = S(200) * 10**-11
    k_11 = S(1500) # Assuming k_22 = k_11

    eps0 = S(9) * 10**-12
    rho = S(4500)

    # For v2 (Longitudinal mode along z)
    e_33 = d_33 * c_33
    eps_33 = k_33 * eps0
    c_prime_33 = c_33 + e_33**2 / eps_33
    v2_sq = c_prime_33 / rho
    
    # For v1 (Transverse mode along z)
    e_24 = d_24 * c_44
    eps_22 = k_11 * eps0
    c_prime_44 = c_44 + e_24**2 / eps_22
    v1_sq = c_prime_44 / rho
    
    # Step 3: Calculate R^2 = (v2/v1)^2
    R_sq = v2_sq / v1_sq

    # The equation for A_max^8 is:
    # (1 - R^4) * X^2 + R^2 * X - beta = 0
    # where X = A_max^8
    
    R_four = R_sq**2
    beta = S('17/324')
    
    # Coefficients of the quadratic equation a*X^2 + b*X + c = 0
    a = 1 - R_four
    b = R_sq
    c = -beta

    # Step 4: Solve the quadratic equation for X = A_max^8
    # Using the quadratic formula: X = (-b + sqrt(b^2 - 4ac)) / 2a
    # We take the positive root because X must be positive.
    discriminant = b**2 - 4*a*c
    # We convert to float for calculation, as sympy's sqrt on large fractions can be slow
    X = (-b + np.sqrt(float(discriminant))) / (2*a)

    # A_max is the maximum amplitude |Φ|
    A_max = float(X)**(S(1)/S(8))

    # The problem asks for (1 - max|Φ|)
    result = 1 - A_max

    # Print the equation with final numerical values
    # We need to print each number in the equation, as requested.
    # Note: A_max^16 = X^2
    print("The algebraic equation for the soliton amplitude A is:")
    print(f"-{beta} + ({R_sq}) * A**8 + ({1 - R_four}) * A**16 = 0")
    print("\nSolving this for A gives max|Φ| = A_max.")
    print(f"The calculated value for A_max is: {A_max}")
    
    print("\nThe final result (1 - max|Φ|) is:")
    print(result)

solve_problem()