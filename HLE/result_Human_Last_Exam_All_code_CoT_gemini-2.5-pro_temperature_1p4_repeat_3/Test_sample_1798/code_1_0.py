import numpy as np
from scipy.integrate import quad

def solve_and_print_results():
    """
    This function solves the problem by finding the circle's equation,
    calculating alpha and the required probability, and printing the results.
    """
    # Step 1: Find the equation of the circle.
    # The circle passes through A=(1,0), B=(10,0), and C=(4,3).
    # The perpendicular bisector of the segment AB is x = (1+10)/2 = 5.5.
    # So, the x-coordinate of the center is h = 5.5.
    # Using the circle equation (x-h)^2 + (y-k)^2 = r^2 with points A and C:
    # For A(1,0): (1-5.5)^2 + (0-k)^2 = r^2  => (-4.5)^2 + k^2 = r^2 => 20.25 + k^2 = r^2
    # For C(4,3): (4-5.5)^2 + (3-k)^2 = r^2  => (-1.5)^2 + (3-k)^2 = r^2 => 2.25 + 9 - 6k + k^2 = r^2
    # Equating the two gives: 20.25 + k^2 = 11.25 - 6k + k^2 => 9 = -6k => k = -1.5
    # Then, r^2 = 20.25 + (-1.5)^2 = 20.25 + 2.25 = 22.5
    h = 5.5
    k = -1.5
    r2 = 22.5

    print("The equation of the circle is (x - h)^2 + (y - k)^2 = r^2.")
    print("The values of the parameters are:")
    print(f"h = {h}")
    print(f"k = {k}")
    print(f"r^2 = {r2}")
    print(f"The equation is: (x - {h})^2 + (y - ({k}))^2 = {r2}")
    
    # Step 2: Define the function f(x)
    # y = k +/- sqrt(r^2 - (x-h)^2). To pass through C(4,3) we must choose the '+' sign.
    # y = -1.5 + sqrt(22.5 - (x-5.5)^2).
    def f(x):
        discriminant = r2 - (x - h)**2
        # Ensure discriminant is non-negative for inputs at the boundary
        discriminant = np.maximum(0, discriminant)
        return k + np.sqrt(discriminant)

    # Step 3: Calculate alpha
    # alpha = 1 / integral of f(x) from 1 to 10.
    total_integral, _ = quad(f, 1, 10)
    alpha = 1.0 / total_integral

    # Step 4: Calculate P(X < 3)
    # P(X < 3) = alpha * integral of f(x) from 1 to 3.
    partial_integral, _ = quad(f, 1, 3)
    probability_X_lt_3 = alpha * partial_integral
    
    print("\nThe calculated values are:")
    print(f"alpha = {alpha}")
    print(f"P(X < 3) = {probability_X_lt_3}")

solve_and_print_results()