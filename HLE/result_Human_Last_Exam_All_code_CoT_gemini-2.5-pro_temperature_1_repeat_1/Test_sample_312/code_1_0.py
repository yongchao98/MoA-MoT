from fractions import Fraction

def solve_frostman_decay():
    """
    Calculates the smallest possible value of c based on Wolff's theorem
    for the decay of Fourier transforms of Frostman measures.
    """
    # The dimension of the Frostman measure is given as s = 8/5.
    s = Fraction(8, 5)
    
    # The problem asks for the smallest possible c such that for every 8/5-Frostman measure mu,
    # the L^2 norm of its Fourier transform on a circle of radius r is O(r^{c+epsilon}).
    # This is a known result in harmonic analysis.
    
    print(f"The dimension of the measure is s = {s.numerator}/{s.denominator}.")
    
    # The critical exponent is given by a theorem from Thomas Wolff, which applies for s > 3/2.
    threshold = Fraction(3, 2)
    print(f"We check the condition s > 3/2: {s} > {threshold} is {s > threshold}.")
    
    # The theorem states that the L^2 norm squared decays with an exponent of -(s - 3/2) / (s - 1).
    print("\nCalculating the exponent for the squared L^2 norm...")
    
    # Numerator of the exponent fraction
    num = s - threshold
    print(f"Numerator of the exponent is: s - 3/2 = {s.numerator}/{s.denominator} - {threshold.numerator}/{threshold.denominator} = {num.numerator}/{num.denominator}")
    
    # Denominator of the exponent fraction
    den = s - 1
    print(f"Denominator of the exponent is: s - 1 = {s.numerator}/{s.denominator} - 1 = {den.numerator}/{den.denominator}")
    
    # Exponent for the squared norm
    squared_norm_exponent = -num / den
    print(f"Exponent for the squared L^2 norm is: -({num}) / ({den}) = {squared_norm_exponent.numerator}/{squared_norm_exponent.denominator}")
    
    # The exponent c is for the L^2 norm, which is the square root of the L^2 norm squared.
    # Therefore, we must divide the exponent by 2.
    c = squared_norm_exponent / 2
    print("\nThe exponent c is half of the exponent for the squared norm.")
    
    # Final equation and result
    print("\nThe final equation for c is:")
    s_num, s_den = s.numerator, s.denominator
    thresh_num, thresh_den = threshold.numerator, threshold.denominator
    c_num, c_den = c.numerator, c.denominator
    
    print(f"c = ( - ( {s_num}/{s_den} - {thresh_num}/{thresh_den} ) / ( {s_num}/{s_den} - 1 ) ) / 2 = {c_num}/{c_den}")

solve_frostman_decay()