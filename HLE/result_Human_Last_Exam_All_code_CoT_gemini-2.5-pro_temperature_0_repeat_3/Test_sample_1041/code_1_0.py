import math

def solve_random_walk_probability():
    """
    Calculates the probability that a 2D simple random walk visits a target set
    before leaving a large disk.
    """
    # Problem parameters
    R = 1000.0
    S_x, S_y = 0.0, 300.0
    A1_x, A1_y = 0.0, 0.0
    A2_x, A2_y = 2.0, 0.0

    # Euler-Mascheroni constant
    gamma = 0.5772156649015329

    # --- Calculations based on the plan ---

    # 1. Determine the center of the target set A
    zA_x = (A1_x + A2_x) / 2.0
    zA_y = (A1_y + A2_y) / 2.0

    # 2. Calculate the capacitary radius of the target set A
    # Distance between the two points in set A
    d = math.sqrt((A2_x - A1_x)**2 + (A2_y - A1_y)**2)
    # Capacitary radius of a single point for SRW on Z^2
    r0 = math.exp(-gamma) / 2.0
    # Effective capacitary radius of the two-point set A
    # r_A = sqrt(d * r_0) = sqrt(2 * e^(-gamma)/2) = e^(-gamma/2)
    rA = math.exp(-gamma / 2.0)

    # 3. Calculate the distance from the starting point S to the center of A
    dist_S_zA = math.sqrt((S_x - zA_x)**2 + (S_y - zA_y)**2)

    # 4. Apply the harmonic measure formula for the hitting probability
    numerator = math.log(R / dist_S_zA)
    denominator = math.log(R / rA)
    probability = numerator / denominator

    # --- Output the results ---

    print("The probability P is calculated using the formula: P ≈ ln(R / |S - z_A|) / ln(R / r_A)")
    print("\nThe final equation with all the initial numbers is:")
    # The formula shows all the components before they are computed.
    print(f"P ≈ ln({R} / sqrt(({S_x} - {zA_x})**2 + ({S_y} - {zA_y})**2)) / ln({R} / exp(-{gamma} / 2.0))")

    print("\nIntermediate calculated values:")
    print(f"|S - z_A| = {dist_S_zA}")
    print(f"r_A = {rA}")

    print("\nSubstituting these values into the formula:")
    print(f"P ≈ ln({R} / {dist_S_zA}) / ln({R} / {rA})")
    print(f"P ≈ {numerator} / {denominator}")
    
    print(f"\nFinal calculated probability: {probability}")
    print(f"\nThe probability rounded to three significant digits is: {probability:.3g}")

solve_random_walk_probability()