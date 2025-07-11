import random
import math

def get_circumcircle(p1, p2, p3):
    """
    Calculates the circumcenter and squared radius of the circle passing through three points.
    Returns None if the points are collinear.
    """
    x1, y1 = p1
    x2, y2 = p2
    x3, y3 = p3

    # Using the formula from Wikipedia for the circumcenter
    D = 2 * (x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2))

    if abs(D) < 1e-10:  # The points are collinear
        return None, None

    # Unpack for cleaner formulas
    x1_sq, y1_sq = x1*x1, y1*y1
    x2_sq, y2_sq = x2*x2, y2*y2
    x3_sq, y3_sq = x3*x3, y3*y3
    
    # Circumcenter coordinates (cx, cy)
    cx = ( (x1_sq + y1_sq) * (y2 - y3) + (x2_sq + y2_sq) * (y3 - y1) + (x3_sq + y3_sq) * (y1 - y2) ) / D
    cy = ( (x1_sq + y1_sq) * (x3 - x2) + (x2_sq + y2_sq) * (x1 - x3) + (x3_sq + y3_sq) * (x2 - x1) ) / D
    
    # Squared radius is the squared distance from the center to any of the points
    radius_sq = (x1 - cx)**2 + (y1 - cy)**2
    
    return (cx, cy), radius_sq

def simulate_duck_problem(num_trials):
    """
    Runs a Monte Carlo simulation for the duck problem.
    """
    inside_count = 0
    valid_trials = 0

    for _ in range(num_trials):
        # Generate three random ducks for the circle
        p1 = (random.uniform(0, 1), random.uniform(0, 1))
        p2 = (random.uniform(0, 1), random.uniform(0, 1))
        p3 = (random.uniform(0, 1), random.uniform(0, 1))
        
        # Generate the fourth duck
        p4 = (random.uniform(0, 1), random.uniform(0, 1))
        
        center, radius_sq = get_circumcircle(p1, p2, p3)
        
        if center is None:
            # Skip this trial if the first three points are collinear
            continue
        
        valid_trials += 1
        
        # Calculate the squared distance from the fourth duck to the circumcenter
        dist_sq = (p4[0] - center[0])**2 + (p4[1] - center[1])**2
        
        # Check if the fourth duck is inside the circle
        if dist_sq < radius_sq:
            inside_count += 1
            
    if valid_trials == 0:
        print("No valid trials were completed (all points were collinear).")
        return

    # Calculate and print the final probability estimate
    probability = inside_count / valid_trials
    
    print("Monte Carlo Simulation Results:")
    print(f"Total Trials: {num_trials}")
    print(f"Valid (non-collinear) Trials: {valid_trials}")
    print(f"Duck #4 was inside the circle: {inside_count} times")
    print("\nFinal Equation:")
    print(f"Probability ≈ {inside_count} / {valid_trials}")
    print(f"Estimated Probability = {probability:.6f}")
    print("\nThe theoretical probability is exactly 0.5.")

if __name__ == "__main__":
    # Number of simulations to run. More trials give a more accurate estimate.
    number_of_trials = 1000000
    simulate_duck_problem(number_of_trials)