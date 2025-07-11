import sympy

def solve_c3_max():
    """
    This function calculates the maximum value of c3 by carrying out the
    maximization of the intermediate function and then computing the final result.
    """
    # In the text above, we derived that we need to maximize the function
    # h(x, y) = 5*x*y*(x-y) for x in (0, 1] and y in (0, 1].
    # Let's use sympy to analyze this function.
    x, y = sympy.symbols('x y', real=True)
    h = 5 * x * y * (x - y)

    # We found that the maximum occurs on the boundary of the domain.
    # Let's verify this. We check the interior for critical points.
    # The gradient is [5*y*(2*x-y), 5*x*(x-2*y)]. Setting it to zero
    # for x,y > 0 gives x=0 and y=0, which is not in the interior.
    # So the maximum must be on the boundary of the square [0,1] x [0,1].

    # Boundary 1: x = 1, for y in [0, 1]
    h_x1 = h.subs(x, 1)  # h(1, y) = 5*y*(1-y)
    # To find the maximum of h(1,y), we take the derivative w.r.t. y and set to 0.
    # d/dy (5y - 5y^2) = 5 - 10y = 0  => y = 1/2.
    y_crit_val = sympy.Rational(1, 2)
    max_val_boundary1 = h_x1.subs(y, y_crit_val)

    # Boundary 2: y = 1, for x in [0, 1]
    h_y1 = h.subs(y, 1) # h(x, 1) = 5*x*(x-1)
    # The maximum value of x^2 - x on [0,1] is 0 (at x=0 and x=1).
    max_val_boundary2 = 0

    # The other boundaries (x=0 or y=0) result in h=0.
    
    # The maximum value of the function h(x,y) is the maximum found on the boundaries.
    max_integral_val = max(max_val_boundary1, max_val_boundary2, 0)
    
    # Now we compute the maximum value of c3.
    # c3_max = (7/2) * max_integral_val
    factor = sympy.Rational(7, 2)
    c3_max = factor * max_integral_val

    print("Step 1: Find the maximum of the intermediate function.")
    print(f"The maximum value of the integral part is: {max_integral_val}")
    print("\nStep 2: Calculate the maximum value of c3.")
    print("The final equation for the maximum value of c3 is:")
    print(f"c3_max = ({factor.p}/{factor.q}) * ({max_integral_val.p}/{max_integral_val.q}) = {c3_max.p}/{c3_max.q}")
    
    print("\nThe numbers in the final equation are:")
    print(f"Numerator 1: {factor.p}")
    print(f"Denominator 1: {factor.q}")
    print(f"Numerator 2: {max_integral_val.p}")
    print(f"Denominator 2: {max_integral_val.q}")
    print(f"Final Numerator: {c3_max.p}")
    print(f"Final Denominator: {c3_max.q}")
    
    print("\nThe maximum value of c3 as a decimal is:")
    print(float(c3_max))

solve_c3_max()