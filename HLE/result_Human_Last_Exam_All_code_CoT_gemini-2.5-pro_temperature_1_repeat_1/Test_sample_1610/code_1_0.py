import sympy

def solve_largest_r():
    """
    Calculates the largest real number r based on the described geometric problem.

    The value r is the ratio of the area of the central region to the
    total integral of a weight function A(x,y) over the 4x4 square.
    """

    # Define the symbolic variable
    x = sympy.symbols('x')

    # Define the 1D weight function A_x(x) = min(3, x) - max(0, x-1) for x in [0, 4]
    # We can define it piecewise for integration.
    # For x in [0, 1]: A_x(x) = x - 0 = x
    # For x in [1, 3]: A_x(x) = x - (x-1) = 1
    # For x in [3, 4]: A_x(x) = 3 - (x-1) = 4 - x

    # Calculate the integral of A_x(x) from 0 to 4
    integral_ax_part1 = sympy.integrate(x, (x, 0, 1))
    integral_ax_part2 = sympy.integrate(1, (x, 1, 3))
    integral_ax_part3 = sympy.integrate(4 - x, (x, 3, 4))

    total_integral_ax = integral_ax_part1 + integral_ax_part2 + integral_ax_part3

    # The total integral of A(x,y) over the 4x4 square is (integral of A_x(x))^2
    total_integral = total_integral_ax ** 2

    # The area of the central region C = [1,3]x[1,3] is (3-1)*(3-1)
    area_C = (3 - 1) * (3 - 1)

    # The largest value for r is the ratio of these two quantities
    r = sympy.Rational(area_C, total_integral)
    
    # Print the result as an equation
    print(f"The largest real number r is the solution to the equation:")
    print(f"r = Area(Central Region) / TotalIntegratedWeight")
    print(f"r = {area_C} / {total_integral} = {r.evalf()}")
    
solve_largest_r()