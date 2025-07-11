import numpy as np
from scipy.integrate import quad

def solve_and_print():
    """
    This function solves the problem by finding the equation of the circle,
    calculating the normalization constant alpha, and then finding the
    probability P(X < 3).
    """

    # Step 1: Define the function f(x) from the circle's equation.
    # The points are A=(1,0), B=(10,0), C=(4,3).
    # The perpendicular bisector of AB is x = (1+10)/2 = 5.5, so the center's x-coordinate h=5.5.
    # By substituting the points into (x-h)^2 + (y-k)^2 = r^2, we find:
    # h = 5.5, k = -1.5, r^2 = 22.5
    # The equation of the circle is (x - 5.5)^2 + (y + 1.5)^2 = 22.5.
    # We solve for y and take the upper branch to pass through C(4,3):
    # y = f(x) = sqrt(22.5 - (x - 5.5)^2) - 1.5
    
    h = 5.5
    k = -1.5
    r_squared = 22.5

    def f(x):
        # Ensure the argument of sqrt is non-negative.
        # This is guaranteed for x in [1,10] as the circle's x-domain is approx [0.75, 10.25].
        return np.sqrt(r_squared - (x - h)**2) + k

    # Step 2: Calculate the total integral to find alpha.
    # alpha = 1 / integral of f(x) from 1 to 10.
    total_integral, _ = quad(f, 1, 10)

    if total_integral <= 0:
        print("Error: The integral of f(x) must be positive.")
        return

    alpha = 1 / total_integral

    # Step 3: Calculate the partial integral for P(X < 3).
    # P(X < 3) = alpha * integral of f(x) from 1 to 3.
    partial_integral, _ = quad(f, 1, 3)

    prob_x_less_than_3 = alpha * partial_integral

    # Output the results
    print("The equation of the circle that passes through A(1,0), B(10,0), and C(4,3) is:")
    print(f"(x - {h})^2 + (y - ({k}))^2 = {r_squared}")
    print(f"(x - 5.5)^2 + (y + 1.5)^2 = 22.5\n")
    
    print("The probability P(X < 3) is calculated as follows:")
    print(f"P(X < 3) = (Integral of f(x) from 1 to 3) / (Integral of f(x) from 1 to 10)")
    
    # As requested, outputting each number in the final equation for the probability
    print(f"P(X < 3) = {partial_integral} / {total_integral}")
    print(f"P(X < 3) = {prob_x_less_than_3}\n")

    print("--- Final Answers ---")
    print(f"The value of alpha is: {alpha}")
    print(f"The value of P(X < 3) is: {prob_x_less_than_3}")
    
    # Storing the final result for the grading format
    global final_answer
    final_answer = f"{alpha}, {prob_x_less_than_3}"


if __name__ == '__main__':
    solve_and_print()
    # The final answer will be printed to the console.
    # To conform to the output format, we would have something like this:
    # print(f"<<<{final_answer}>>>")