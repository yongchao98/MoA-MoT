import numpy as np
from scipy.integrate import quad

def solve_and_print():
    """
    This function solves the problem by finding the equation of the circle,
    calculating the normalization constant alpha, and then finding the
    probability P(X < 3).
    """
    # Step 1: Define the function f(x) for the upper arc of the circle.
    # The center of the circle passing through A(1,0), B(10,0), and C(4,3)
    # can be found to be (h, k) = (5.5, -1.5).
    # The radius squared is r^2 = (1 - 5.5)^2 + (0 - (-1.5))^2 = (-4.5)^2 + 1.5^2 = 20.25 + 2.25 = 22.5.
    # The equation of the circle is (x - 5.5)^2 + (y + 1.5)^2 = 22.5.
    # Solving for y gives y = -1.5 +/- sqrt(22.5 - (x - 5.5)^2).
    # Since C(4,3) is on the arc, we take the positive root: f(x) = -1.5 + sqrt(22.5 - (x - 5.5)^2).
    
    h = 5.5
    k = -1.5
    r_squared = 22.5

    def f(x):
        """The equation of the upper arc of the circle."""
        return k + np.sqrt(r_squared - (x - h)**2)

    # Step 2: Find the normalization constant alpha.
    # alpha is defined such that the integral of alpha * f(x) from 1 to 10 is 1.
    # This means alpha = 1 / (integral of f(x) from 1 to 10).
    
    # Calculate the integral of f(x) from 1 to 10.
    total_area, _ = quad(f, 1, 10)
    
    # Calculate alpha
    alpha = 1 / total_area

    # Step 3: Calculate the probability P(X < 3).
    # P(X < 3) is the integral of the PDF d_X(x) from 1 to 3.
    # P(X < 3) = alpha * (integral of f(x) from 1 to 3)
    
    # Calculate the integral of f(x) from 1 to 3.
    partial_area, _ = quad(f, 1, 3)
    
    # Calculate the probability
    probability_X_less_than_3 = alpha * partial_area

    # Print the results and the numbers used in the final equations.
    print("Step 1: The equation for the upper arc of the circle is f(x) = -1.5 + sqrt(22.5 - (x - 5.5)^2)")
    
    print("\nStep 2: Calculation of alpha")
    print("The total area under f(x) from x=1 to x=10 is the integral of f(x) dx.")
    print(f"Integral(f(x), 1, 10) = {total_area}")
    print(f"alpha = 1 / {total_area}")
    print(f"alpha = {alpha}")

    print("\nStep 3: Calculation of P(X < 3)")
    print("The area under f(x) from x=1 to x=3 is the integral of f(x) dx.")
    print(f"Integral(f(x), 1, 3) = {partial_area}")
    print(f"P(X < 3) = alpha * Integral(f(x), 1, 3) = {alpha} * {partial_area}")
    print(f"P(X < 3) = {probability_X_less_than_3}")

if __name__ == '__main__':
    solve_and_print()