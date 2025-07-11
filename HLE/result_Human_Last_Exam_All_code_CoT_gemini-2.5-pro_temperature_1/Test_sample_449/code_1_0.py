import numpy as np

def solve_random_walk_probability():
    """
    Calculates the probability that a 2D SRW conditioned to avoid the origin
    also avoids the origin's four neighbors.
    """
    # 1. Define constants and the starting point
    x0_coords = (3000, 4000)
    norm_x0 = np.sqrt(x0_coords[0]**2 + x0_coords[1]**2)
    
    # 2. Calculate the potential kernel at a neighbor of the origin, y in B.
    # a(1,0) = 4/pi - 1. This is the same for all 4 neighbors by symmetry.
    a_B = 4 / np.pi - 1
    
    # 3. Calculate the potential kernel at the starting point x0 using the asymptotic formula.
    # a(x) ~ (2/pi) * (log|x| + gamma + log(sqrt(8)))
    gamma = np.euler_gamma  # Euler-Mascheroni constant
    a_x0 = (2 / np.pi) * (np.log(norm_x0) + gamma + np.log(np.sqrt(8)))
    
    # 4. Calculate the probability of the conditioned walk HITting the neighbors.
    # This is approximated by a_B / a_x0.
    prob_hit = a_B / a_x0
    
    # 5. The probability of NEVER hitting the neighbors is 1 minus the hitting probability.
    prob_never_hit = 1 - prob_hit
    
    # 6. Print the results and the equation
    print("This script calculates the probability for a conditioned 2D random walk.")
    print(f"The starting point is x0={x0_coords}, with distance from origin |x0| = {norm_x0:.0f}.")
    print("The probability of the walk *never* hitting the origin's neighbors is P = 1 - P_hit.")
    print("P_hit is the probability of hitting the neighbors, approximated by a(y) / a(x0), where y is a neighbor.")
    
    print("\n--- Calculation ---")
    print(f"Potential kernel at a neighbor, a(y) = 4/pi - 1 = {a_B:.5f}")
    print(f"Potential kernel at the start point, a(x0) = {a_x0:.5f}")
    
    print("\nFinal Equation:")
    print(f"P = 1 - (a(y) / a(x0))")
    print(f"P = 1 - {a_B:.5f} / {a_x0:.5f}")
    print(f"P = 1 - {prob_hit:.5f}")
    print(f"P = {prob_never_hit:.5f}")
    
    print(f"\nThe approximate answer with two significant digits is: {prob_never_hit:.2f}")

solve_random_walk_probability()
<<<0.96>>>