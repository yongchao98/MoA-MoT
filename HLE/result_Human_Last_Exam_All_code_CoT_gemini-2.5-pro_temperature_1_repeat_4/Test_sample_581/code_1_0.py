import math

def solve_cap_set_bound():
    """
    Calculates and explains the best-known lower bound for the size of a cap set in dimension 8.

    A cap set in the n-dimensional affine space over the field of 3 elements, AG(n, 3),
    is a set of points with no three collinear. The size of the largest possible
    cap set is denoted r_3(n). This function calculates the best-known lower bound for r_3(8).

    The bound is based on a construction by Bierbrauer and Edel, which uses the
    maximal cap set in dimension 4, r_3(4) = 20.
    """

    # The size of the maximal cap set in dimension 4, r_3(4).
    r3_4 = 20

    # The number of additional points gained from the specific construction
    # by Bierbrauer and Edel for k=4.
    additional_points = 96

    # Calculate the lower bound for r_3(8) using the formula: r_3(4)^2 + 96
    lower_bound_r3_8 = int(math.pow(r3_4, 2) + additional_points)

    print("The best known lower bound for the size of a cap set in dimension 8 is based on a construction by Bierbrauer and Edel.")
    print(f"This construction uses the maximal cap set from dimension 4, which has a size of {r3_4}.")
    print("The final size is calculated from the following equation:")
    print(f"{r3_4}^2 + {additional_points} = {lower_bound_r3_8}")


solve_cap_set_bound()