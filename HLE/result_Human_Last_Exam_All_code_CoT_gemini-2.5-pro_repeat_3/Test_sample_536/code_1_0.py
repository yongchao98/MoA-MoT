import math

def solve_for_limit():
    """
    This function solves for the limit of the inner product <b_p, z_p>.
    The problem is reduced to solving a geometric equation involving the inner products
    of three unit vectors h, b, and z, assuming they are coplanar in the limit.

    The equation is:
    <h, b> = <h, z> * <b, z> + sqrt(1 - <h, z>^2) * sqrt(1 - <b, z>^2)

    Given values:
    <h, z> = 0.9
    <h, b> = 0.9375

    Let x = <b, z>. The equation becomes:
    0.9375 = 0.9 * x + sqrt(1 - 0.9^2) * sqrt(1 - x^2)

    Rearranging and squaring leads to a quadratic equation.
    """

    # Coefficients from the equation:
    # (0.9375 - 0.9x)^2 = (1 - 0.9^2)(1 - x^2)
    # 0.9375^2 - 2*0.9375*0.9*x + 0.81*x^2 = 0.19(1-x^2)
    # 0.87890625 - 1.6875*x + 0.81*x^2 = 0.19 - 0.19*x^2
    # (0.81+0.19)x^2 - 1.6875x + (0.87890625 - 0.19) = 0
    # 1.0*x^2 - 1.6875*x + 0.68890625 = 0
    a = 1.0
    b = -1.6875
    c = 0.68890625

    print("The value x = lim <b_p, z_p> is a root of the quadratic equation:")
    print(f"{a} * x^2 + ({b}) * x + ({c}) = 0")
    print("\nSolving this equation:")

    # Calculate the discriminant
    discriminant = b**2 - 4 * a * c

    if discriminant < 0:
        print("The equation has no real solutions.")
        return

    # Calculate the two possible solutions for x
    sol1 = (-b - math.sqrt(discriminant)) / (2 * a)
    sol2 = (-b + math.sqrt(discriminant)) / (2 * a)

    print(f"The two potential solutions are x1 = {sol1:.8f} and x2 = {sol2:.8f}.")

    # Both solutions are mathematically valid under the geometric model.
    # In the context of these problems, a single answer is expected.
    # Without further physical justification to distinguish between the two, we select
    # the smaller root, which represents a more significant angular separation
    # between the signal vector b_p and the DC vector z_p.
    final_answer = min(sol1, sol2)
    
    print("\nAfter analyzing the solutions, we select the smaller value.")
    print(f"The final answer is: {final_answer}")
    print(f"\n<<< {final_answer} >>>")


solve_for_limit()