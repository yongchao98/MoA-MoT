import math

def solve():
    """
    Calculates the probability that a 2D simple random walk, conditioned to avoid the origin,
    starting from (0,1), hits the neighbors of (3600,0).
    """
    
    # --- Step 1: Define constants ---
    pi = math.pi
    gamma = math.gamma
    
    # --- Step 2: Calculate the potential kernel at the starting point x0 = (0,1) ---
    # a(0,1) is exactly known to be 4/pi
    a_x0 = 4 / pi
    
    # --- Step 3: Calculate the potential kernel at the target point z = (3600,0) ---
    # We use the asymptotic formula for a(z) for large |z|.
    # a(z) ~ (2/pi) * ln(|z|) + (2*gamma + ln(8))/pi
    
    # Parameters for the formula
    z_norm = 3600.0
    log_z = math.log(z_norm)
    log_8 = math.log(8.0)
    
    # The constant term in the asymptotic expansion
    K = (2 * gamma + log_8) / pi
    
    # Calculate a(z)
    a_z = (2 / pi) * log_z + K
    
    # --- Step 4: Calculate the final probability using the refined formula ---
    # P(x0) approx a(x0) / (2*a(z) - a(x0))
    
    numerator = a_x0
    denominator = 2 * a_z - a_x0
    
    probability = numerator / denominator
    
    # --- Step 5: Print the results ---
    print("This problem asks for the probability P that a conditioned random walk starting from x0=(0,1) hits the set A (neighbors of z=(3600,0)).")
    print("This is approximately P(x0) = a(x0) / (2*a(z) - a(x0)), where a(x) is the potential kernel.")
    print("\nCalculation steps:")
    print(f"1. Potential at starting point: a(0,1) = 4/pi = {a_x0:.4f}")
    print(f"2. Potential at target point: a(3600,0) approx (2/pi)*ln(3600) + K = {a_z:.4f}")
    print(f"   (where K = (2*gamma + ln(8))/pi = {K:.4f})")
    
    print("\nFinal Equation:")
    print(f"P = {numerator:.4f} / (2 * {a_z:.4f} - {a_x0:.4f})")
    print(f"P = {numerator:.4f} / ({denominator:.4f})")
    print(f"P = {probability:.4f}")
    
    print(f"\nThe approximate probability with two significant digits is: {probability:.2g}")

solve()
<<<0.11>>>