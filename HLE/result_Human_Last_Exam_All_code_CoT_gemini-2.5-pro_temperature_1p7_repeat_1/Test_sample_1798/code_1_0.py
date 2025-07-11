import numpy as np
from scipy.integrate import quad

def solve_and_print_results():
    """
    Solves the problem by finding the equation of the circle,
    calculating the normalization constant alpha, and computing the probability P(X < 3).
    """
    # Step 1: Define the function f(x) based on the circle's properties.
    # The center (h,k) and radius squared (r_sq) are derived from the points A(1,0), B(10,0), C(4,3).
    # h = (1+10)/2 = 5.5
    # (1-5.5)^2 + (0-k)^2 = r^2  => 20.25 + k^2 = r^2
    # (4-5.5)^2 + (3-k)^2 = r^2  => 2.25 + (3-k)^2 = r^2
    # 20.25 + k^2 = 2.25 + 9 - 6k + k^2 => 9 = -6k => k = -1.5
    # r^2 = 20.25 + (-1.5)^2 = 22.5
    h = 5.5
    k = -1.5
    r_sq = 22.5

    def f(x):
        """
        Represents the y-coordinate of the upper semi-circle.
        f(x) = k + sqrt(r^2 - (x-h)^2)
        """
        return k + np.sqrt(r_sq - (x - h)**2)

    # Step 2: Calculate the normalization constant alpha.
    # alpha = 1 / integral of f(x) from 1 to 10.
    total_area, _ = quad(f, 1, 10)
    alpha = 1 / total_area

    # Step 3: Calculate the probability P(X < 3).
    # P(X < 3) = alpha * integral of f(x) from 1 to 3.
    area_lt_3, _ = quad(f, 1, 3)
    prob_lt_3 = alpha * area_lt_3
    
    # Step 4: Output the results as requested.
    print("--- Problem Analysis ---")
    print("The circle passes through A=(1,0), B=(10,0), and C=(4,3).")
    print("The probability density function is d_X(x) = alpha * f(x) on [1, 10], where y=f(x) is the upper arc of the circle.")
    print("\n--- Equation of the Circle ---")
    print(f"Center (h,k) = ({h}, {k})")
    print(f"Radius squared r^2 = {r_sq}")
    # Printing each number in the equation as requested
    print(f"The equation is: (x - {h})^2 + (y - ({k}))^2 = {r_sq}")
    
    print("\n--- Calculation of alpha ---")
    print(f"The integral of f(x) from 1 to 10 is: {total_area:.8f}")
    print(f"alpha = 1 / {total_area:.8f}")
    print(f"The value of alpha is: {alpha:.8f}")
    
    print("\n--- Calculation of P(X < 3) ---")
    print(f"The integral of f(x) from 1 to 3 is: {area_lt_3:.8f}")
    print(f"P(X < 3) = alpha * (integral from 1 to 3)")
    print(f"The value of P(X < 3) is: {prob_lt_3:.8f}")
    
    # Final answer for direct extraction if needed
    # print(f"\nFinal values: alpha = {alpha}, P(X < 3) = {prob_lt_3}")

solve_and_print_results()