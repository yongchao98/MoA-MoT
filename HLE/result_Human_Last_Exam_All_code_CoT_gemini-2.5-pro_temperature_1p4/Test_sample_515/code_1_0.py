import math

def solve_probability():
    """
    Calculates the probability that a 2D random walk starting at (0,1),
    conditioned to never visit the origin, hits the neighbors of (3600,0).

    This is equivalent to finding the probability that a standard simple random
    walk starting at (0,1) hits the neighbors of (3600,0) before returning to the origin.
    The probability p is given by the ratio of the potential kernel a(x) at the
    starting point x and the target point z: p ~ a(x) / a(z).
    """

    # Start and target points
    x_coord = (0, 1)
    z_coord = (3600, 0)
    
    # Use literature value for a(1,0) = 4/pi
    # a(0,1) is the same by symmetry.
    a_x = 4 / math.pi
    
    # Use the asymptotic formula for a(3600,0)
    # a(z) ~ (2/pi)*ln(||z||) + (2*gamma + ln(8))/pi
    # where ||z|| is the distance from the origin.
    z_norm = 3600
    gamma = 0.5772156649  # Euler-Mascheroni constant

    # Numerator of the expression for the probability
    numerator = 4
    
    # Denominator
    denominator = 2 * math.log(z_norm) + 2 * gamma + math.log(8)
    
    # Calculate the probability
    probability = numerator / denominator
    
    # Print the explanation and the final equation
    print("The problem asks for the probability that a 2D random walk, started at x=(0,1) and conditioned to avoid the origin, reaches the neighbors of z=(3600,0).")
    print("This probability 'p' can be approximated by the ratio of the potential kernel 'a' evaluated at the start and target locations:")
    print("p ≈ a(x) / a(z)")
    print("\nUsing known formulas and values for the potential kernel:")
    print("a(x) = a(0,1), which by symmetry is a(1,0). A known value from literature is a(1,0) = 4/π.")
    print("a(z) = a(3600,0), which can be found using the asymptotic formula:")
    print("a(z) ≈ (2/π) * ln(||z||) + (2γ + ln(8))/π")
    print("where ||z|| = 3600 and γ is the Euler-Mascheroni constant.")
    print("\nThe ratio simplifies to:")
    print("p ≈ (4/π) / [ (2/π)ln(3600) + (2γ + ln(8))/π ]")
    print("p ≈ 4 / [ 2*ln(3600) + 2γ + ln(8) ]")
    print("\nPlugging in the values:")
    print(f"p ≈ {numerator} / [ 2*ln({z_norm}) + 2*{gamma:.4f} + ln(8) ]")
    print(f"p ≈ {numerator} / [ {2 * math.log(z_norm):.4f} + {2 * gamma:.4f} + {math.log(8):.4f} ]")
    print(f"p ≈ {numerator} / [ {denominator:.4f} ]")
    print(f"p ≈ {probability:.4f}")
    
    # Final answer with two significant digits
    print(f"\nThe approximate probability with two significant digits is: {probability:.2f}")

solve_probability()
<<<0.20>>>