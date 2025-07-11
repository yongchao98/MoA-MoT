import math

def solve_rw_probability():
    """
    Calculates the probability for a 2D random walk to visit a set
    before leaving a large disk.
    """
    # Parameters from the problem statement
    z = (0, 300)
    a1 = (0, 0)
    a2 = (2, 0)
    R = 1000
    r0 = 1.0  # Effective radius of a single lattice point

    # Calculate the required distances
    dist_z_a1 = math.sqrt((z[0] - a1[0])**2 + (z[1] - a1[1])**2)
    dist_z_a2 = math.sqrt((z[0] - a2[0])**2 + (z[1] - a2[1])**2)
    dist_a1_a2 = math.sqrt((a1[0] - a2[0])**2 + (a1[1] - a2[1])**2)

    # Calculate the terms in the formula using natural logarithm (ln)
    log_term1_num = math.log(R / dist_z_a1)
    log_term2_num = math.log(R / dist_z_a2)
    
    log_term1_den = math.log(R / r0)
    log_term2_den = math.log(R / dist_a1_a2)

    # Calculate numerator and denominator
    numerator = log_term1_num + log_term2_num
    denominator = log_term1_den + log_term2_den

    # Calculate the final probability
    probability = numerator / denominator

    # Print the breakdown of the calculation as requested
    print("The probability P(z) is estimated using the formula:")
    print("P(z) approx = [ln(R/|z-a1|) + ln(R/|z-a2|)] / [ln(R/r0) + ln(R/|a1-a2|)]")
    print("\nPlugging in the values:")
    print(f"z = {z}, a1 = {a1}, a2 = {a2}, R = {R}, r0 = {r0:.1f}")
    print(f"|z-a1| = {dist_z_a1:.3f}")
    print(f"|z-a2| = {dist_z_a2:.3f}")
    print(f"|a1-a2| = {dist_a1_a2:.3f}")
    
    print("\nCalculating the final equation:")
    print(f"P approx = [ln({R}/{dist_z_a1:.3f}) + ln({R}/{dist_z_a2:.3f})] / [ln({R}/{r0:.1f}) + ln({R}/{dist_a1_a2:.1f})]")
    print(f"P approx = [{log_term1_num:.4f} + {log_term2_num:.4f}] / [{log_term1_den:.4f} + {log_term2_den:.4f}]")
    print(f"P approx = {numerator:.4f} / {denominator:.4f}")
    
    print(f"\nThe calculated probability is: {probability:.5f}")
    print(f"The probability rounded to three significant digits is: {probability:.3g}")

solve_rw_probability()