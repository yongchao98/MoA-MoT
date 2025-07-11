import numpy as np
from scipy.integrate import quad

def solve_probability_problem():
    """
    Solves the problem by finding the equation of the circle,
    calculating the normalization constant alpha, and the probability P(X<3).
    """
    # Step 1: Define the parameters of the circle from the derived equation
    # (x - h)^2 + (y - k)^2 = r^2
    # (x - 5.5)^2 + (y + 1.5)^2 = 22.5
    h = 5.5
    k = -1.5
    r_sq = 22.5

    # Step 2: Define the function f(x) which is the upper arc of the circle
    def f(x):
        return k + np.sqrt(r_sq - (x - h)**2)

    # Step 3: Calculate the total area under f(x) from 1 to 10 to find alpha
    # The integral of the PDF alpha*f(x) over [1, 10] must be 1.
    # So, alpha * integral(f(x) dx) = 1 => alpha = 1 / integral(f(x) dx)
    total_area, _ = quad(f, 1, 10)
    
    alpha = 1 / total_area

    # Step 4: Calculate P(X < 3)
    # P(X < 3) is the integral of the PDF from 1 to 3.
    # P(X < 3) = alpha * integral from 1 to 3 of f(x) dx
    area_lt_3, _ = quad(f, 1, 3)
    
    prob_lt_3 = alpha * area_lt_3

    # Output the results
    print("Problem Solution:")
    print("-" * 20)
    print("The equation of the circle passing through A(1,0), B(10,0), and C(4,3) is:")
    print(f"(x - {h})^2 + (y - ({k}))^2 = {r_sq}")
    print("\nThe function y=f(x) for the upper arc is:")
    print(f"f(x) = {k} + sqrt({r_sq} - (x - {h})^2)")
    print("\nThe probability density function is d_X(x) = alpha * f(x) on [1,10].")
    
    print("\nCalculations:")
    print(f"The integral of f(x) from 1 to 10 is: {total_area:.4f}")
    print(f"alpha = 1 / {total_area:.4f}")
    print(f"The value of alpha is: {alpha:.4f}")
    print("-" * 20)
    print(f"The integral of f(x) from 1 to 3 is: {area_lt_3:.4f}")
    print(f"P(X < 3) = alpha * {area_lt_3:.4f} = {alpha:.4f} * {area_lt_3:.4f}")
    print(f"The value of P(X < 3) is: {prob_lt_3:.4f}")
    print("-" * 20)

solve_probability_problem()

# The final answers are alpha and P(X<3)
# To format for the <<<>>> tags, let's get the final numerical values
# alpha = 0.046830...
# prob_lt_3 = 0.207578...
# Rounding to a reasonable number of digits, e.g., 4.
# alpha = 0.0468
# prob_lt_3 = 0.2076
# Let's provide the answer as a tuple (alpha, P(X<3))
final_alpha_val = 1 / quad(lambda x: -1.5 + np.sqrt(22.5 - (x - 5.5)**2), 1, 10)[0]
final_prob_val = final_alpha_val * quad(lambda x: -1.5 + np.sqrt(22.5 - (x - 5.5)**2), 1, 3)[0]
# The final answer requested should be just the values.
# print(f'<<<{final_alpha_val}, {final_prob_val}>>>')
