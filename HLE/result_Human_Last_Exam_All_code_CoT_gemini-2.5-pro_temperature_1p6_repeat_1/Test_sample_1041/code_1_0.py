import math

def solve_random_walk_probability():
    """
    Calculates the probability that a 2D simple random walk visits a target set
    before leaving a disk.
    """
    # Parameters of the problem
    R = 1000.0
    z0 = (0, 300)
    z1 = (0, 0)
    z2 = (2, 0)

    # Calculate the required distances
    d01 = math.sqrt((z0[0] - z1[0])**2 + (z0[1] - z1[1])**2)
    d02 = math.sqrt((z0[0] - z2[0])**2 + (z0[1] - z2[1])**2)
    d12 = math.sqrt((z1[0] - z2[0])**2 + (z1[1] - z2[1])**2)
    d11 = 1.0  # Effective distance for self-potential on a lattice

    # Approximate the Green's functions using ln(R/d)
    # The proportionality constant cancels out in the ratio.
    G_z0_z1 = math.log(R / d01)
    G_z0_z2 = math.log(R / d02)
    G_z1_z1 = math.log(R / d11)
    G_z1_z2 = math.log(R / d12)

    # The probability is the ratio of potentials
    numerator = G_z0_z1 + G_z0_z2
    denominator = G_z1_z1 + G_z1_z2 # By symmetry G_z2_z2=G_z1_z1, G_z2_z1=G_z1_z2
    
    probability = numerator / denominator

    # Output the equation and the result
    print("The probability is approximated by the formula:")
    print("P = [ln(R/|z0-z1|) + ln(R/|z0-z2|)] / [ln(R/|z1-z1|) + ln(R/|z1-z2|)]")
    print("\nPlugging in the values:")
    print(f"R = {R}")
    print(f"|z0-z1| = {d01:.4f}")
    print(f"|z0-z2| = {d02:.4f}")
    print(f"|z1-z1| = {d11:.4f} (lattice spacing)")
    print(f"|z1-z2| = {d12:.4f}")
    
    print(f"\nP = [ln({R}/{d01:.4f}) + ln({R}/{d02:.4f})] / [ln({R}/{d11:.4f}) + ln({R}/{d12:.4f})]")
    print(f"P = [{G_z0_z1:.4f} + {G_z0_z2:.4f}] / [{G_z1_z1:.4f} + {G_z1_z2:.4f}]")
    print(f"P = {numerator:.4f} / {denominator:.4f}")
    
    print(f"\nFinal Probability: {probability:.3f}")

solve_random_walk_probability()