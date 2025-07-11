import math

def solve():
    """
    This function solves the geometry problem step-by-step as outlined above.
    """
    # From the geometric analysis, we found several key parameters.
    # The ratio k0 = BE/EC
    k0 = 1/2
    # The length of the equal sides AB and AC, denoted by L
    L = 3
    
    # We found the cosine of angle BAC using the Law of Cosines on triangle AGH
    # AG = 6, AH = 2, GH = 5
    AG = 6
    AH = 2
    GH = 5
    # GH^2 = AG^2 + AH^2 - 2 * AG * AH * cos(BAC)
    # 25 = 36 + 4 - 2 * 6 * 2 * cos(BAC)
    # 25 = 40 - 24 * cos(BAC)
    # 15 = 24 * cos(BAC)
    cos_BAC = 15/24

    # From Stewart's Theorem on triangle ABC with cevian AE, we got L^2 = AE^2 + BE*CE
    # We need to find BE and CE.
    # BC^2 = L^2 + L^2 - 2*L*L*cos(BAC)
    BC_sq = 2 * L**2 * (1 - cos_BAC)
    # BC_sq = 2 * 9 * (1 - 15/24) = 18 * (9/24) = 18 * 3/8 = 54/8 = 27/4
    # So BC = sqrt(27)/2 = 3*sqrt(3)/2
    
    # BE/EC = k0 = 1/2, so BC = BE + EC = BE + 2*BE = 3*BE
    # BE = BC/3
    # CE = 2*BE
    BC = math.sqrt(27/4)
    BE = BC / 3
    CE = 2 * BE
    BE_CE = BE * CE
    
    # Now find AE^2 from L^2 = AE^2 + BE*CE
    AE_sq = L**2 - BE_CE
    
    # AE^2 = 9 - ( (3*sqrt(3)/6) * (2*3*sqrt(3)/6) )
    # BE = 3*sqrt(3)/6 = sqrt(3)/2
    # CE = 3*sqrt(3)/3 = sqrt(3)
    # BE*CE = (sqrt(3)/2)*sqrt(3) = 3/2
    calculated_BE_CE = 3/2
    AE_sq_calc = 9 - calculated_BE_CE
    
    # AE_sq = 15/2. So AE = sqrt(15/2) = sqrt(30)/2
    
    a = 1
    b = 30
    c = 2
    
    print(f"The length of AE is sqrt({AE_sq_calc.numerator*2}/{AE_sq_calc.denominator*2}) = sqrt({b})/{c}")
    print(f"This can be written as (a*sqrt(b))/c with a = {a}, b = {b}, c = {c}.")

    result = a + b + c
    print(f"The value of a+b+c is {a} + {b} + {c} = {result}.")

solve()