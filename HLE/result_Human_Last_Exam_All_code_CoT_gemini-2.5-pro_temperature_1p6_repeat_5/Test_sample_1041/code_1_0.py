import math

def solve_random_walk_probability():
    """
    Calculates the probability that a 2D simple random walk visits a target set
    before leaving a large disk.
    """
    # Parameters from the problem statement
    # Starting point is S = (0, 300)
    # Target set is A = {(0,0), (2,0)}
    # Escape disk has radius R = 1000, centered at the origin.

    R = 1000.0  # Radius of the escape disk
    S = 300.0   # Distance from origin to the starting point
    d = 2.0     # Distance between the two points in the target set A

    # The probability P is approximated by the formula:
    # P = (2 * ln(R/S)) / (2 * ln(R) - ln(d))
    
    # Calculate numerator and denominator
    numerator = 2 * math.log(R / S)
    denominator = 2 * math.log(R) - math.log(d)
    
    # Calculate the probability
    probability = numerator / denominator

    # Print the explanation and the equation with its values
    print("The problem is to find the probability that a 2D random walk starting at (0, 300)")
    print("visits the set {(0,0), (2,0)} before leaving the disk of radius 1000.")
    print("\nThis probability can be approximated using potential theory for random walks.")
    print("The final formula for the probability P is:")
    print("P = (2 * ln(R / S)) / (2 * ln(R) - ln(d))")
    print("\nWith the following values:")
    print(f"R (escape radius) = {R}")
    print(f"S (starting distance) = {S}")
    print(f"d (target distance) = {d}")

    # Print the equation with the numbers plugged in
    print("\nPlugging in the numbers, the equation is:")
    print(f"P = (2 * ln({R} / {S})) / (2 * ln({R}) - ln({d}))")

    # Print the final result
    print(f"\nThe calculated probability is: {probability}")
    print(f"\nThe probability rounded to three significant digits is: {probability:.3g}")

solve_random_walk_probability()
<<<0.183>>>