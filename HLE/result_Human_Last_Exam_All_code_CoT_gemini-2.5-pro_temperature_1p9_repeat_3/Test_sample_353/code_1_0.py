import sympy
import numpy as np

def solve_bdf4_stability_angle():
    """
    Finds the exact value of the A(alpha)-stability angle for the BDF4 method.
    """
    # Step 1: Solve the polynomial for x = cos(theta)
    # The condition for the tangent from the origin to the stability boundary
    # of BDF4 reduces to a polynomial equation in x = cos(theta).
    x = sympy.Symbol('x')
    poly_eq = 5*x**4 - 16*x**3 + 18*x**2 - 8*x + 1
    
    # We find the roots of this polynomial.
    # We can use numerical solvers to find the roots first, then verify them symbolically.
    p_coeffs = [5, -16, 18, -8, 1]
    roots = np.roots(p_coeffs)
    
    # The roots are approximately [1., 1., 1., 0.2].
    # Let's verify the rational root x=1/5.
    cos_theta_val = sympy.Rational(1, 5)

    print("Step 1: The polynomial equation for x = cos(theta) is:")
    print(f"   {poly_eq} = 0")
    print(f"The roots are 1 (a triple root) and 1/5.")
    print(f"The root corresponding to the maximal angle is x = {cos_theta_val}.")
    print("-" * 20)
    
    # Step 2: Calculate sin(theta)
    # We take the positive root, as the region is symmetric with respect to the real axis.
    sin_theta_val = sympy.sqrt(1 - cos_theta_val**2)
    print(f"Step 2: From cos(theta) = {cos_theta_val}, we get sin(theta) = {sin_theta_val}.")
    print("-" * 20)

    # Step 3: Calculate the coordinates (X, Y) of the point z(theta)
    # The real and imaginary parts of z(theta) depend on cos(k*theta) and sin(k*theta).
    # Using De Moivre's formula:
    # cos(2t) = 2cos^2(t)-1, sin(2t) = 2sin(t)cos(t)
    # cos(3t) = 4cos^3(t)-3cos(t), sin(3t) = 3sin(t)-4sin^3(t) = sin(t)(3-4sin^2(t))=sin(t)(4cos^2(t)-1)
    # cos(4t) = 8cos^4(t)-8cos^2(t)+1, sin(4t) = 2sin(2t)cos(2t)
    
    c1 = cos_theta_val
    s1 = sin_theta_val
    c2 = 2*c1**2 - 1
    s2 = 2*s1*c1
    c3 = 4*c1**3 - 3*c1
    s3 = s1*(4*c1**2-1)
    c4 = 2*c2**2 - 1 # More stable than using c1^4
    s4 = 2*s2*c2
    
    # BDF4 coefficients for z(theta)
    a0 = sympy.Rational(25, 12)
    a1 = -4
    a2 = 3
    a3 = -sympy.Rational(4, 3)
    a4 = sympy.Rational(1, 4)

    # X = Re(z), Y = Im(z)
    X = a0*1 + a1*c1 + a2*c2 + a3*c3 + a4*c4
    Y = - (a1*s1 + a2*s2 + a3*s3 + a4*s4)
    
    print("Step 3: Calculating the coordinates (X, Y) of the tangent point z.")
    print(f"   X = {X}")
    print(f"   Y = {Y}")
    print("-" * 20)

    # Step 4: Calculate tan(alpha) = -Y/X
    tan_alpha = -Y/X

    print("Step 4: The stability angle alpha is given by arctan(-Y/X).")
    print(f"   tan(alpha) = -({Y}) / ({X}) = {tan_alpha}")
    
    # Final simplified result for tan(alpha)
    num, den = tan_alpha.as_numer_denom()
    
    # We want to format the result as `arctan(A*sqrt(B)/C)`
    sqrt_term = 1
    if isinstance(num, sympy.Mul):
      for term in num.args:
        if isinstance(term, sympy.sqrt):
          sqrt_term = term.args[0]
          break
    num_coeff = sympy.simplify(num / sympy.sqrt(sqrt_term))

    print("-" * 20)
    print("Final Result:")
    print("The exact stability angle alpha for BDF4 is:")
    print(f"alpha = arctan({num_coeff} * sqrt({sqrt_term}) / {den})")

solve_bdf4_stability_angle()