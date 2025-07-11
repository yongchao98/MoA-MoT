import math

def solve_random_walk_problem():
    """
    Calculates the probability that a 2D simple random walk visits a target set
    before leaving a large disk, using a potential theory approximation.
    """

    # Parameters from the problem statement
    R = 1000
    z = (0, 300)
    w1 = (0, 0)
    w2 = (2, 0)

    # Helper function for Euclidean distance
    def distance(p1, p2):
        return math.sqrt((p1[0] - p2[0])**2 + (p1[1] - p2[1])**2)

    # Calculate the required distances
    d_z_w1 = distance(z, w1)
    d_z_w2 = distance(z, w2)
    d_w1_w2 = distance(w1, w2)
    # The effective distance for a point to itself (self-interaction) is
    # approximated by the lattice spacing, which is 1.
    d_w1_w1 = 1.0

    # Calculate the logarithmic terms for the formula
    log_R_dzw1 = math.log(R / d_z_w1)
    log_R_dzw2 = math.log(R / d_z_w2)
    log_R_dw1w1 = math.log(R / d_w1_w1)
    log_R_dw1w2 = math.log(R / d_w1_w2)

    # Calculate the numerator and denominator of the probability formula
    numerator = log_R_dzw1 + log_R_dzw2
    denominator = log_R_dw1w1 + log_R_dw1w2

    # Final probability
    prob = numerator / denominator

    # Output the steps of the calculation
    print("The probability P is approximated by the formula derived from potential theory:")
    print("P = (log(R/|z-w1|) + log(R/|z-w2|)) / (log(R/|w1-w1|) + log(R/|w1-w2|))")
    print("\nSubstituting the given values:")
    print(f"R = {R}")
    print(f"z = {z}, w1 = {w1}, w2 = {w2}")
    print(f"|z-w1| = {d_z_w1:.3f}")
    print(f"|z-w2| = {d_z_w2:.3f}")
    print(f"|w1-w1| is approximated by the lattice spacing = {d_w1_w1}")
    print(f"|w1-w2| = {d_w1_w2:.3f}")
    
    print("\nThe equation becomes:")
    print(f"P = (log({R}/{d_z_w1:.3f}) + log({R}/{d_z_w2:.3f})) / (log({R}/{d_w1_w1:.1f}) + log({R}/{d_w1_w2:.1f}))")
    print(f"P = (log({R/d_z_w1:.4f}) + log({R/d_z_w2:.4f})) / (log({R/d_w1_w1:.1f}) + log({R/d_w1_w2:.1f}))")
    print(f"P = ({log_R_dzw1:.5f} + {log_R_dzw2:.5f}) / ({log_R_dw1w1:.5f} + {log_R_dw1w2:.5f})")
    print(f"P = {numerator:.5f} / {denominator:.5f}")
    
    print("\nThe final calculated probability is:")
    print(f"P = {prob:.5f}")

solve_random_walk_problem()