import math

def calculate_cross_section_ratio():
    """
    Calculates the ratio of differential cross-sections for scattering 
    off a magnetic monopole vs. a magnetic dipole.
    """
    
    # The derived formula for the ratio R is:
    # R = (3 * e_m^2 * hbar^2) / (p^2 * mu^2 * sin(theta)^2)
    # After substituting the given values:
    # e_m = e/16
    # mu = 25 * mu_B = 25 * (e*hbar / (2*m*c))
    # p = m*v = m * (c/100)
    # The physical constants (e, hbar, m, c) cancel out, yielding:
    # R = 3 / (4 * sin(theta)^2), where theta = pi/30
    # Wait, my previous derivation was R = (3 * e_m^2 * hbar^2) / (p^2 * mu^2 * sin^2(theta))
    # Let me re-check my math from the thought process.
    # R = (3 * e_m^2 * hbar^2) / (4 * p^2 * mu^2 * sin^2(theta/2) * cos^2(theta/2))
    # Using sin(theta) = 2*sin(theta/2)cos(theta/2), we get 4*sin^2(theta/2)cos^2(theta/2) = sin^2(theta)
    # So R = (3 * e_m^2 * hbar^2) / (p^2 * mu^2 * sin^2(theta)) is correct.
    # Substituting the values for e_m, mu, p:
    # Numerator term: 3 * (e^2/256) * hbar^2
    # Denominator term: (m^2*c^2/10000) * (625 * e^2*hbar^2 / (4*m^2*c^2)) * sin^2(theta)
    # Denominator simplifies: (625 / 40000) * e^2 * hbar^2 * sin^2(theta) = (1/64) * e^2 * hbar^2 * sin^2(theta)
    # Ratio R = (3 * e^2*hbar^2 / 256) / ( (1/64) * e^2 * hbar^2 * sin^2(theta) )
    # R = (3/256) / ( (1/64) * sin^2(theta) ) = (3/256) * (64 / sin^2(theta)) = 3 / (4 * sin^2(theta))
    # The formula is correct.

    theta = math.pi / 30
    
    numerator = 3.0
    denominator_factor = 4.0
    
    # Print the equation with the numbers plugged in
    print(f"The final equation for the ratio is:")
    print(f"Ratio = {numerator} / ({denominator_factor} * sin^2(pi / 30))")
    
    # Calculate the value
    ratio_value = numerator / (denominator_factor * math.sin(theta)**2)
    
    # Print the calculated result
    print("\nThe calculated ratio is:")
    print(ratio_value)
    
    # Output the final answer in the required format
    print(f"\n<<<{ratio_value}>>>")

calculate_cross_section_ratio()