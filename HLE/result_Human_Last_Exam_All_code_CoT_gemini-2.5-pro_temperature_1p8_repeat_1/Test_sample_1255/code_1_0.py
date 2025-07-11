import math

def check_nodal_curve_parametrization(t):
    """
    This function demonstrates the parametrization of the smooth part of the
    nodal cubic curve y^2 = x^2 * (x + 1).

    The smooth points of this curve form a group isomorphic to the multiplicative group.
    A parameter 't' from the multiplicative group maps to a point (x, y) on the curve.
    The parametrization used here is:
    x = t^2 - 1
    y = t * (t^2 - 1)

    Args:
        t (float or int): A non-zero value representing an element in the multiplicative group.
    """
    if t == 0:
        print("The parameter 't' cannot be zero.")
        return

    # Calculate (x, y) coordinates from the parameter t
    x = t**2 - 1
    y = t * (t**2 - 1)

    # The equation of the nodal cubic curve is y^2 = x^2 * (x + 1)
    lhs = y**2
    rhs = x**2 * (x + 1)

    print(f"For parameter t = {t}:")
    print(f"The corresponding point on the curve is (x, y) = ({x}, {y}).")
    print(f"Now, let's verify this point satisfies the curve equation: y^2 = x^2 * (x + 1).")
    print(f"Left Hand Side (LHS) = y^2 = ({y})^2 = {lhs}")
    print(f"Right Hand Side (RHS) = x^2 * (x + 1) = ({x})^2 * (({x}) + 1) = {rhs}")

    # Check if LHS and RHS are close enough for floating point numbers
    if math.isclose(lhs, rhs):
        print("Verification successful: The point lies on the curve.\n")
    else:
        print("Verification failed: The point does not lie on the curve.\n")

# Demonstrate with a few values for t
check_nodal_curve_parametrization(2)
check_nodal_curve_parametrization(3)
check_nodal_curve_parametrization(-0.5)
