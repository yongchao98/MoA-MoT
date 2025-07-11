import math

def estimate_mass():
    """
    Calculates and prints the steps for estimating the rock's mass
    using fractions with integers no larger than 7.
    """
    # Problem values approximated with fractions using integers <= 7
    # Approximate density (0.9) with 6/7
    density_num, density_den = 6, 7
    # Use standard 4/3 from the volume formula
    vol_const_num, vol_const_den = 4, 3
    # Approximate pi with 3/1
    pi_num, pi_den = 3, 1
    # Original radius (0.5) is 1/2
    radius_num, radius_den = 1, 2

    print("Yes, the parrot can estimate the mass.")
    print("Here are the instructions for the calculation:")
    print("\nMass = (density) * (4/3) * (pi) * (radius)^3")
    print(f"We will use the following approximations:")
    print(f"density ≈ {density_num}/{density_den}")
    print(f"pi ≈ {pi_num}/{pi_den}")
    print(f"radius = {radius_num}/{radius_den}")
    
    # Perform the calculation
    # Numerator calculation
    final_num = density_num * vol_const_num * pi_num * (radius_num**3)
    # Denominator calculation
    final_den = density_den * vol_const_den * pi_den * (radius_den**3)
    
    print("\nThe calculation is:")
    print(f"Mass ≈ ({density_num}/{density_den}) * ({vol_const_num}/{vol_const_den}) * ({pi_num}/{pi_den}) * ({radius_num}/{radius_den})³")
    print(f"Mass ≈ ({density_num} * {vol_const_num} * {pi_num} * {radius_num**3}) / ({density_den} * {vol_const_den} * {pi_den} * {radius_den**3})")
    print(f"Mass ≈ {final_num} / {final_den}")

    # Simplify the final fraction
    common_divisor = math.gcd(final_num, final_den)
    simplified_num = final_num // common_divisor
    simplified_den = final_den // common_divisor
    
    print(f"\nSimplified, the estimated mass is {simplified_num}/{simplified_den} kg.")

estimate_mass()