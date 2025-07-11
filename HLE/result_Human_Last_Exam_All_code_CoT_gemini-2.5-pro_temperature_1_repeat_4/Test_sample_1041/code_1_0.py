import math

def solve_random_walk_probability():
    """
    Calculates the probability that a 2D random walk visits a target set
    before leaving a large disk.
    """
    # Problem parameters
    R = 1000.0
    z0 = (0.0, 300.0)
    a1 = (0.0, 0.0)
    a2 = (2.0, 0.0)

    # --- Step 1: Calculate the terms for the numerator ---
    # The numerator is log(R/|z0-a1|) + log(R/|z0-a2|)
    d_z0_a1 = math.hypot(z0[0] - a1[0], z0[1] - a1[1])
    d_z0_a2 = math.hypot(z0[0] - a2[0], z0[1] - a2[1])
    
    term1_num = math.log(R / d_z0_a1)
    term2_num = math.log(R / d_z0_a2)
    numerator = term1_num + term2_num

    # --- Step 2: Calculate the terms for the denominator ---
    # The denominator is log(R/r_e) + log(R/|a1-a2|)
    
    # Calculate effective radius r_e
    gamma = 0.5772156649015328  # Euler-Mascheroni constant
    neg_log_re = (gamma + math.log(8)) / 2.0
    
    term1_den = math.log(R) + neg_log_re

    # Calculate distance between a1 and a2
    d_a1_a2 = math.hypot(a1[0] - a2[0], a1[1] - a2[1])
    term2_den = math.log(R / d_a1_a2)
    
    denominator = term1_den + term2_den

    # --- Step 3: Calculate the final probability ---
    probability = numerator / denominator
    
    # --- Step 4: Print the results ---
    print("This script calculates the probability that a 2D random walk starting at (0, 300)")
    print("visits the set {(0,0), (2,0)} before leaving the disk of radius 1000.")
    print("\nThe calculation is based on the formula: P = num / den")
    print("where num = log(R/|z0-a1|) + log(R/|z0-a2|)")
    print("and   den = log(R/r_e) + log(R/|a1-a2|)\n")

    print("--- Calculated Values ---")
    print(f"log(R/|z0-a1|) = log({R:.0f}/{d_z0_a1:.3f}) = {term1_num:.4f}")
    print(f"log(R/|z0-a2|) = log({R:.0f}/{d_z0_a2:.3f}) = {term2_num:.4f}")
    print(f"Numerator = {term1_num:.4f} + {term2_num:.4f} = {numerator:.4f}\n")
    
    print(f"log(R/r_e) = {term1_den:.4f}")
    print(f"log(R/|a1-a2|) = log({R:.0f}/{d_a1_a2:.3f}) = {term2_den:.4f}")
    print(f"Denominator = {term1_den:.4f} + {term2_den:.4f} = {denominator:.4f}\n")
    
    print("--- Final Equation ---")
    print(f"P = {numerator:.4f} / {denominator:.4f}\n")

    print("--- Result ---")
    print(f"The probability is approximately: {probability:.5f}")
    print(f"The probability to three significant digits is: {probability:.3g}")


if __name__ == '__main__':
    solve_random_walk_probability()