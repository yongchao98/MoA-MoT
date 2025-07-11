import math

def display_transformation_x():
    """
    This function prints the general finite transformation for the x-coordinate,
    which leaves the given PDE invariant. The transformation group depends on
    arbitrary constants C_T, C_A, and C_B derived from the Lie symmetry analysis.
    """
    
    # The general transformation is derived by integrating the infinitesimal generators.
    # The infinitesimals tau and xi are found to be:
    #   tau = C_T
    #   xi = C_A + C_B * exp(k1*t)
    # where C_T, C_A, C_B are arbitrary constants.
    
    # Integrating d(x_bar)/d(epsilon) = xi(t_bar) where d(t_bar)/d(epsilon) = tau(t_bar)
    # yields two cases for the finite transformation x_bar(epsilon).

    print("The general representation for the Lie group transformation on the spatial coordinate x,")
    print("denoted as x_bar, depends on the parameter for time translation (C_T).\n")
    print("-" * 80)

    # Case 1: The time translation parameter C_T is non-zero.
    print("Case 1: For a non-zero time translation (C_T != 0):")
    case1_eq = "x_bar(t, x; epsilon) = x + C_A * epsilon + (C_B * exp(k1 * t) / (k1 * C_T)) * (exp(k1 * C_T * epsilon) - 1)"
    print(case1_eq)
    print("\nThe number appearing in this equation is -1.")
    print("-" * 80)
    
    # Case 2: The time translation parameter C_T is zero.
    # This case is the limit of Case 1 as C_T -> 0.
    print("Case 2: For a zero time translation (C_T = 0):")
    case2_eq = "x_bar(t, x; epsilon) = x + (C_A + C_B * exp(k1 * t)) * epsilon"
    print(case2_eq)
    print("\nThere are no explicit numbers in this equation besides the implied '1's as coefficients.")

    print("-" * 80)
    print("Where:")
    print("  - x_bar is the transformed spatial coordinate.")
    print("  - (t, x) are the original coordinates.")
    print("  - epsilon is the parameter of the Lie group transformation.")
    print("  - k1 is the constant from the heat equation's source term.")
    print("  - C_T, C_A, C_B are arbitrary constants from the linear combination of the infinitesimal generators.")

if __name__ == '__main__':
    display_transformation_x()
