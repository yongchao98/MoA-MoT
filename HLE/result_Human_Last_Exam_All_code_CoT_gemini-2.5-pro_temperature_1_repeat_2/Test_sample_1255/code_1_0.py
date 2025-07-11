import sympy

def check_singularity():
    """
    This function demonstrates the singularity of a curve, which is relevant to
    the counterexample for the main question.

    The question asks if the underlying scheme of a log group G is always a group scheme.
    The answer is No. A key counterexample is the logarithmic multiplicative group
    over a log point S. Its underlying scheme G_sch is not a group scheme because
    its special fiber, given by the equation xy = 0, is not a group scheme.

    A connected group scheme over a field must be smooth (non-singular).
    The curve xy=0 is singular at the origin (0,0), which prevents it from
    being a group scheme.

    This script verifies the singularity at (0,0) by checking the partial derivatives
    of the defining polynomial F(x, y) = xy. A point is singular if all partial
    derivatives are zero at that point.
    """

    # Define the variables
    x, y = sympy.symbols('x y')

    # Define the equation of the special fiber
    F = x * y

    # Calculate the partial derivatives
    F_dx = sympy.diff(F, x)
    F_dy = sympy.diff(F, y)

    # Point to check for singularity
    point = {x: 0, y: 0}

    # Evaluate the partial derivatives at the point
    F_dx_at_point = F_dx.subs(point)
    F_dy_at_point = F_dy.subs(point)

    print(f"The equation for the special fiber is F(x, y) = {F} = 0.")
    print(f"The point to check for singularity is (x, y) = ({point[x]}, {point[y]}).")
    print(f"The partial derivative with respect to x is dF/dx = {F_dx}.")
    print(f"The partial derivative with respect to y is dF/dy = {F_dy}.")
    print(f"Evaluating dF/dx at the point: {F_dx_at_point}")
    print(f"Evaluating dF/dy at the point: {F_dy_at_point}")

    if F_dx_at_point == 0 and F_dy_at_point == 0:
        print("\nConclusion: Both partial derivatives are 0 at the point (0,0).")
        print("Therefore, the point is a singular point.")
        print("Since the special fiber is a singular connected scheme, it cannot be a group scheme.")
        print("This shows that the underlying scheme of the log multiplicative group is not a group scheme.")
    else:
        print("\nThe point is not singular.")

# Run the check
check_singularity()
