import math

def solve_and_compute_coeff():
    """
    This function solves for the optimal cost per bit for Bob's problem.
    """
    # We need to find the positive real root of the equation r^3 - r^2 - 1 = 0.
    # We can use Newton's method for a quick and accurate approximation.
    # f(r) = r^3 - r^2 - 1
    # f'(r) = 3*r^2 - 2*r
    # An initial guess, knowing the root is between 1 and 2.
    r = 1.5 
    # Iterating 10 times gives high precision.
    for _ in range(10):
        r = r - (r**3 - r**2 - 1) / (3*r**2 - 2*r)

    # The minimal cost per bit K is given by 1 / log_2(r).
    # This is calculated as ln(2) / ln(r).
    cost_coefficient = math.log(2) / math.log(r)
    
    print("The minimal cost is expressed as C = K * log2(n!), where K is the cost per bit.")
    print("K is determined by the root of the equation derived from the costs of Type 2 questions (1 for 'yes', 3 for 'no').")
    print("The equation for the base 'r' is: 1*r^3 - 1*r^2 + 0*r - 1 = 0")
    print(f"The positive real root of this equation is r ≈ {r:.6f}")
    print("The final equation for the minimal cost per bit K is: K = 1 / log2(r)")
    print(f"The calculated value for K is: {cost_coefficient:.3f}")

solve_and_compute_coeff()
