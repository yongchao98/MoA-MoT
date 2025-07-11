import numpy as np
from scipy.integrate import quad

def solve_circle_pdf():
    """
    Solves for the normalization constant alpha and a probability for a PDF
    derived from a circle passing through three points.
    """
    # Step 1: Find the equation of the circle y = f(x).
    # The three points are A=(1,0), B=(10,0), and C=(4,3).
    # The perpendicular bisector of the segment AB gives the x-coordinate of the center.
    h = (1.0 + 10.0) / 2.0
    
    # Using the general equation of a circle (x-h)^2 + (y-k)^2 = r^2 and the points,
    # we can set up a system of equations.
    # (1-h)^2 + (0-k)^2 = r^2   => (1-5.5)^2 + k^2 = r^2   => 20.25 + k^2 = r^2
    # (4-h)^2 + (3-k)^2 = r^2   => (4-5.5)^2 + (3-k)^2 = r^2 => 2.25 + (3-k)^2 = r^2
    #
    # Equating the two expressions for r^2:
    # 20.25 + k^2 = 2.25 + 9 - 6k + k^2
    # 20.25 = 11.25 - 6k
    # 9 = -6k
    k = -1.5
    
    # Substitute k back to find r^2
    r_sq = 20.25 + k**2
    
    # The equation of the circle is (x - 5.5)^2 + (y + 1.5)^2 = 22.5
    # We solve for y to get y = f(x). We choose the positive square root because
    # point C=(4,3) has a positive y-value relative to the center's y-coordinate.
    # y = -1.5 + sqrt(22.5 - (x-5.5)^2)
    
    def f(x):
        """Equation of the upper part of the circle passing through the points."""
        return k + np.sqrt(r_sq - (x - h)**2)

    # Step 2: Find alpha
    # The PDF d_X(x) = alpha * f(x) must integrate to 1 over its domain [1, 10].
    # alpha * integral from 1 to 10 of f(x) dx = 1
    # We calculate the integral of f(x) from 1 to 10 numerically.
    total_integral, _ = quad(f, 1, 10)
    
    # alpha is the reciprocal of this integral.
    alpha = 1.0 / total_integral
    
    # Step 3: Calculate P(X < 3)
    # P(X < 3) is the integral of the PDF from the domain's start (1) to 3.
    # P(X < 3) = alpha * integral from 1 to 3 of f(x) dx
    integral_1_to_3, _ = quad(f, 1, 3)
    
    prob_X_lt_3 = alpha * integral_1_to_3

    # Output the results
    print("The equation of the circle that passes through A, B, and C is:")
    print(f"(x - {h})^2 + (y - ({k}))^2 = {r_sq}")
    print("\nThe function y = f(x) derived from the circle is:")
    print(f"y = {k} + sqrt({r_sq} - (x - {h})^2)")
    print("\n--- Results ---")
    print(f"The value of alpha is: {alpha}")
    print(f"The value of P(X < 3) is: {prob_X_lt_3}")

if __name__ == '__main__':
    solve_circle_pdf()