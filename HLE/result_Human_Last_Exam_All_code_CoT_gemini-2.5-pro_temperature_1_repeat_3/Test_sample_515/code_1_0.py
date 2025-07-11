import math

def calculate_probability():
    """
    This function calculates the desired probability based on the derived formula.
    p(x_0) approx (a(x_0) + a(z) - a(x_0 - z)) / (2 * a(x_0))
    where x_0 = (0,1) and z = (3600,0).
    """

    x0 = (0, 1)
    z = (3600, 0)
    x0_minus_z = (x0[0] - z[0], x0[1] - z[1])
    
    # The value of the potential kernel a(x) for x=(0,1) or (1,0) is known to be 1/pi.
    a_x0 = 1.0 / math.pi
    
    # For large |x|, the potential kernel a(x) can be approximated by:
    # a(x) ~ (2/pi) * log(|x|) + C
    # We only need the difference a(z) - a(x0-z), so the constant C cancels.
    # a(z) - a(x0-z) = a(3600,0) - a(-3600,1) = a(3600,0) - a(3600,1) by symmetry.
    
    norm_z_sq = z[0]**2 + z[1]**2
    norm_x0_minus_z_sq = x0_minus_z[0]**2 + x0_minus_z[1]**2
    
    # Using the log approximation for the difference
    # a(z) - a(x0-z) approx (1/pi) * (log(norm_z_sq) - log(norm_x0_minus_z_sq))
    # This is -(1/pi) * log(norm_x0_minus_z_sq / norm_z_sq)
    # which is -(1/pi) * log((3600^2 + 1) / 3600^2) = -(1/pi) * log(1 + 1/3600^2)
    
    log_term = math.log(1 + 1.0 / (3600**2))
    diff_a = -log_term / math.pi
    
    # The formula for the probability is:
    # p = (a(x0) + diff_a) / (2 * a(x0))
    
    p = (a_x0 + diff_a) / (2 * a_x0)
    
    # Let's print the components of the final equation
    # The equation is p = (a(0,1) + a(3600,0) - a(3600,1)) / (2 * a(0,1))
    # We approximate this as p = (1 - log(1 + 1/3600^2)) / 2
    
    print("The final probability is calculated using the formula:")
    print("p = (a(0,1) + a(3600,0) - a(3600,1)) / (2 * a(0,1))")
    print("\nThis can be simplified using known values and approximations:")
    print("p ≈ (1 - log(1 + 1/3600^2)) / 2")
    
    term1 = 1.0
    term2 = math.log(1.0 + 1.0/3600**2)
    denominator = 2.0
    
    print(f"\nBreaking down the calculation p ≈ ({term1} - {term2}) / {denominator}")
    
    final_prob = (term1 - term2) / denominator
    
    print(f"\nThe calculated probability is: {final_prob}")
    print(f"The approximate answer with two significant digits is: {final_prob:.2f}")

calculate_probability()