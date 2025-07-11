import math

def solve_random_walk_problem():
    """
    Calculates the probability that a 2D simple random walk starting at (0, 300)
    visits the set {(0,0), (2,0)} before leaving the disk of radius 1000.
    """

    # Problem parameters
    z0 = (0, 300)
    p1 = (0, 0)
    p2 = (2, 0)
    R = 1000.0

    # The effective radius of a point target.
    # We assume the standard continuum approximation where this is the lattice spacing.
    r0 = 1.0

    # Step 1: Calculate the necessary distances
    d1 = math.sqrt((z0[0] - p1[0])**2 + (z0[1] - p1[1])**2)
    d2 = math.sqrt((z0[0] - p2[0])**2 + (z0[1] - p2[1])**2)
    d12 = math.sqrt((p1[0] - p2[0])**2 + (p1[1] - p2[1])**2)

    # Step 2: Calculate the terms for the probability formula based on potential theory.
    # The formula for the probability P is: (log(R/d1) + log(R/d2)) / (log(R/r0) + log(R/d12))
    log_R_d1 = math.log(R / d1)
    log_R_d2 = math.log(R / d2)
    log_R_r0 = math.log(R / r0)
    log_R_d12 = math.log(R / d12)
    
    # Step 3: Calculate the numerator and denominator of the probability formula.
    numerator = log_R_d1 + log_R_d2
    denominator = log_R_r0 + log_R_d12

    # Step 4: Compute the final probability.
    prob = numerator / denominator

    # Output the detailed calculation steps.
    print(f"The starting point is z0 = {z0}")
    print(f"The target set is A = {{{p1}, {p2}}}")
    print(f"The radius of the disk is R = {R}")
    print(f"The effective radius of a point target is assumed to be r0 = {r0}")
    
    print("\n--- Calculation Steps ---")
    
    print("\nStep 1: Distances calculation")
    print(f"Distance from start z0 to p1: d1 = |({z0[0]},{z0[1]}) - ({p1[0]},{p1[1]})| = {d1:.4f}")
    print(f"Distance from start z0 to p2: d2 = |({z0[0]},{z0[1]}) - ({p2[0]},{p2[1]})| = {d2:.4f}")
    print(f"Distance between targets p1 and p2: d12 = |({p1[0]},{p1[1]}) - ({p2[0]},{p2[1]})| = {d12:.4f}")
    
    print("\nStep 2: Logarithmic potential terms")
    print("Formula: P = (log(R/d1) + log(R/d2)) / (log(R/r0) + log(R/d12))")
    print(f"log(R/d1) = log({R}/{d1:.4f}) = {log_R_d1:.4f}")
    print(f"log(R/d2) = log({R}/{d2:.4f}) = {log_R_d2:.4f}")
    print(f"log(R/r0) = log({R}/{r0}) = {log_R_r0:.4f}")
    print(f"log(R/d12) = log({R}/{d12:.4f}) = {log_R_d12:.4f}")
    
    print("\nStep 3: Numerator and Denominator of P")
    print(f"Numerator = {log_R_d1:.4f} + {log_R_d2:.4f} = {numerator:.4f}")
    print(f"Denominator = {log_R_r0:.4f} + {log_R_d12:.4f} = {denominator:.4f}")

    print("\nStep 4: Final Probability Calculation")
    print(f"P = {numerator:.4f} / {denominator:.4f} = {prob:.5f}")

    print(f"\nThe final probability rounded to three significant digits is: {prob:.3g}")

solve_random_walk_problem()