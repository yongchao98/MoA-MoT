import math
from scipy.integrate import quad

def solve_and_print_results():
    """
    Solves the problem by finding the circle equation, the PDF, alpha, and the required probability.
    """
    # Step 1 & 2: Determine the equation of the circle and the function f(x).
    # The points are A=(1,0), B=(10,0), C=(4,3).
    # Solving the system of equations for the circle (x-h)^2 + (y-k)^2 = r^2 yields:
    # (1-h)^2 + k^2 = r^2
    # (10-h)^2 + k^2 = r^2
    # (4-h)^2 + (3-k)^2 = r^2
    # From the first two equations: 1 - 2h = 100 - 20h => 18h = 99 => h = 5.5
    # Substituting h=5.5 into the equations gives k = -1.5 and r^2 = 22.5.
    h = 5.5
    k = -1.5
    r_squared = 22.5
    
    # The circle equation is (x - 5.5)^2 + (y + 1.5)^2 = 22.5.
    # We solve for y to get f(x). Point C(4,3) determines the sign.
    # 3 = -1.5 +/- sqrt(22.5 - (4-5.5)^2) = -1.5 +/- 4.5.
    # 3 = -1.5 + 4.5, so we take the positive root.
    # f(x) = -1.5 + sqrt(22.5 - (x-5.5)^2)
    def f(x):
        radicand = r_squared - (x - h)**2
        # This check is for safety, but limits [1,10] ensure radicand >= 0.
        if radicand < 0:
            return 0.0
        return k + math.sqrt(radicand)

    print("Step 1: Find the equation of the circle that passes through A=(1,0), B=(10,0), and C=(4,3).")
    print(f"The equation is (x - {h})^2 + (y - ({k}))^2 = {r_squared}")
    print(f"(x - 5.5)^2 + (y + 1.5)^2 = 22.5\n")

    # Step 3 & 4: Find the normalization constant alpha.
    # We need to calculate the integral of f(x) from 1 to 10.
    # The integral represents the area under the curve y=f(x) between x=1 and x=10.
    total_area, _ = quad(f, 1, 10)
    
    # alpha is the reciprocal of this area.
    alpha = 1.0 / total_area
    
    print("Step 2: Find the value of alpha.")
    print(f"The probability density function is d_X(x) = alpha * f(x) for x in [1, 10].")
    print("For d_X(x) to be a valid PDF, the integral over the domain must be 1.")
    print(f"Integral of f(x) from 1 to 10 = {total_area}")
    print(f"alpha = 1 / {total_area} = {alpha}\n")

    # Step 5: Calculate P(X < 3).
    # This is the integral of the PDF from 1 to 3.
    # The domain of X is [1, 10], so P(X < 3) is the same as P(1 <= X < 3).
    partial_area, _ = quad(f, 1, 3)
    
    probability = alpha * partial_area
    
    print("Step 3: Find the value of P(X < 3).")
    print(f"P(X < 3) = alpha * (Integral of f(x) from 1 to 3)")
    print(f"Integral of f(x) from 1 to 3 = {partial_area}")
    print(f"P(X < 3) = {alpha} * {partial_area} = {probability}\n")
    
    print("Final Answer:")
    print(f"alpha = {alpha}")
    print(f"P(X < 3) = {probability}")

solve_and_print_results()