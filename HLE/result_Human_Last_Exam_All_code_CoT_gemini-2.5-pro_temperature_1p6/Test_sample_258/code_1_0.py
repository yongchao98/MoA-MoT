import math

def solve():
    """
    Calculates the minimal and maximal numbers of grid's cells that a
    circumference of radius 500 can cross, under the given conditions.
    """
    R = 500

    # According to results from geometry of numbers, for a circle with an integer radius R,
    # the maximum number of crossed cells occurs when the center is at a half-integer
    # coordinate, e.g., (0.5, 0.5). The number of cells is 8R.
    max_cells = 8 * R

    # The minimum number of cells occurs when the center is infinitesimally close to an
    # integer coordinate, e.g., (epsilon, epsilon) for a very small epsilon.
    # The number of cells is 8R - 4.
    min_cells = 8 * R - 4
    
    # We can verify that these configurations satisfy the problem's constraints:
    # 1. Non-tangency: Since R is an integer, a center at (k+0.5, l+0.5) or (k+eps, l+eps)
    #    ensures x_c and y_c are not integers, so the circle is not tangent.
    # 2. No crossing at intersections:
    #    For center (0.5, 0.5), we need (k-0.5)^2 + (l-0.5)^2 != 500^2.
    #    (2k-1)^2/4 + (2l-1)^2/4 != 500^2 => (2k-1)^2 + (2l-1)^2 != 1000^2.
    #    The LHS is a sum of two odd squares, which has the form 4m+2.
    #    The RHS is divisible by 4. So they can never be equal.
    #    For a center near an integer point, it's always possible to choose it
    #    to avoid passing through any grid intersection.
    
    print(f"The radius of the circle is R = {R}.")
    print(f"The maximal number of cells crossed is 8 * R = 8 * {R} = {max_cells}.")
    print(f"The minimal number of cells crossed is 8 * R - 4 = 8 * {R} - 4 = {min_cells}.")
    print(f"The answer is ({min_cells}, {max_cells}).")

solve()
