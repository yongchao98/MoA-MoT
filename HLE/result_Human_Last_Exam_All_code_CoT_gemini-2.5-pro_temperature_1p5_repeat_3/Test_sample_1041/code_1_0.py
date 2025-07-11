import math

def solve_random_walk_probability():
    """
    Calculates the probability for a 2D random walk to hit a target set
    before escaping a large disk, using a potential theory approximation.
    """
    # Problem parameters
    R = 1000.0  # Radius of the outer boundary disk
    S = (0, 300) # Starting point of the walk
    a1 = (0, 0)  # First target point
    a2 = (2, 0)  # Second target point
    
    # In the continuous approximation of a discrete random walk, the effective
    # radius of a single point (to handle the log singularity) is taken
    # as the lattice spacing.
    epsilon = 1.0 
    
    # Distance between the two target points
    d = math.sqrt((a1[0] - a2[0])**2 + (a1[1] - a2[1])**2)

    # Distances from the starting point S to the target points a1 and a2
    dist_S_a1 = math.sqrt((S[0] - a1[0])**2 + (S[1] - a1[1])**2)
    dist_S_a2 = math.sqrt((S[0] - a2[0])**2 + (S[1] - a2[1])**2)

    # The probability P(z) is given by P(z) = mu1*log(R/|z-a1|) + mu2*log(R/|z-a2|).
    # We solve for mu1 and mu2 using the boundary conditions P(a1)=1 and P(a2)=1.
    # 1 = mu1*log(R/epsilon) + mu2*log(R/d)
    # 1 = mu1*log(R/d)     + mu2*log(R/epsilon)
    # By symmetry, mu1 = mu2. Let's call the common value mu.
    # 1 = mu * (log(R/epsilon) + log(R/d))
    # So, the probability is P(S) = mu * (log(R/|S-a1|) + log(R/|S-a2|))
    
    # Numerator of the final probability expression
    numerator_val = math.log(R / dist_S_a1) + math.log(R / dist_S_a2)
    
    # Denominator of the final probability expression
    denominator_val = math.log(R / epsilon) + math.log(R / d)
    
    # Final probability
    probability = numerator_val / denominator_val

    # Print the breakdown of the calculation as requested
    print("The probability P is calculated using the following formula derived from potential theory:")
    print(f"P = [log(R/|S-a1|) + log(R/|S-a2|)] / [log(R/epsilon) + log(R/d)]\n")
    print("Substituting the given values:")
    print(f"R = {R}")
    print(f"S = {S}, a1 = {a1}, a2 = {a2}")
    print(f"|S-a1| = {dist_S_a1:.4f}")
    print(f"|S-a2| = {dist_S_a2:.4f}")
    print(f"epsilon = {epsilon}")
    print(f"d = |a1-a2| = {d}\n")
    
    print("The final equation is:")
    print(f"P = [log({R}/{dist_S_a1:.4f}) + log({R}/{dist_S_a2:.4f})] / [log({R}/{epsilon}) + log({R}/{d})]")
    print(f"P = [{math.log(R / dist_S_a1):.4f} + {math.log(R / dist_S_a2):.4f}] / [{math.log(R / epsilon):.4f} + {math.log(R / d):.4f}]")
    print(f"P = {numerator_val:.4f} / {denominator_val:.4f}")
    print(f"P = {probability:.5f}\n")
    
    print(f"The probability rounded to three significant digits is: {probability:.3g}")

solve_random_walk_probability()