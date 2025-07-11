import math

def solve():
    """
    Calculates the probability that a 2D simple random walk, conditioned on never entering the origin,
    starting from (0,1), will eventually hit the set of the four neighbours of (3600,0).
    """

    # --- Step 1: Define constants ---
    gamma = 0.5772156649
    pi = math.pi
    log2 = math.log(2)
    
    start_point_norm = 1.0
    target_point_norm = 3600.0

    # --- Step 2: Use the asymptotic formula for the potential kernel a(x) ---
    # a(x) approx (2/pi) * log(||x||) + C for large ||x||
    # C = (2*gamma + log(8)) / pi = (2*gamma + 3*log(2)) / pi
    C = (2 * gamma + 3 * log2) / pi

    # --- Step 3: Calculate a(z) for z = (3600, 0) ---
    # ||z|| = 3600, which is large enough for the asymptotic formula.
    log_target_norm = math.log(target_point_norm)
    a_target = (2 / pi) * log_target_norm + C
    
    # --- Step 4: Calculate a(x_0) for x_0 = (0, 1) ---
    # ||x_0|| = 1 is not large. We use the harmonic property of a(x).
    # a(0,1) = (1/4) * [a(0,0) + a(1,1) + a(-1,1) + a(0,2)]
    # a(0,0) = 0 by definition. By symmetry, a(1,1) = a(-1,1).
    # a(0,1) = (1/4) * [2*a(1,1) + a(0,2)]
    # We estimate a(1,1) and a(0,2) using the asymptotic formula, as ||(1,1)|| = sqrt(2) and ||(0,2)|| = 2 are reasonably large.
    # a(1,1) where ||(1,1)|| = sqrt(2)
    a_1_1 = (2 / pi) * math.log(math.sqrt(2)) + C
    # a_1_1 = (1 / pi) * log(2) + C
    
    # a(0,2) where ||(0,2)|| = 2
    a_0_2 = (2 / pi) * math.log(2) + C
    
    # Now, combine them to find a(0,1)
    a_start = (1/4) * (2 * a_1_1 + a_0_2)
    # This simplifies to:
    # a_start = (1/4) * (2 * ((1/pi)*log2 + C) + (2/pi)*log2 + C)
    # a_start = (1/4) * ((4/pi)*log2 + 3*C)
    # a_start = (1/pi)*log2 + (3/4)*C

    # --- Step 5: Compute the probability p ≈ a(x_0)/a(z) ---
    prob = a_start / a_target

    # --- Step 6: Print the results ---
    print(f"The potential kernel constant C is approximately: {C:.4f}")
    print(f"The potential kernel value at the starting point (0,1), a(0,1), is approximately: {a_start:.4f}")
    print(f"The potential kernel value at the target point (3600,0), a(3600,0), is approximately: {a_target:.4f}")
    print(f"The probability is the ratio of these two values:")
    print(f"P = a(0,1) / a(3600,0)")
    print(f"P = {a_start:.4f} / {a_target:.4f}")
    print(f"P ≈ {prob:.2g}")

solve()