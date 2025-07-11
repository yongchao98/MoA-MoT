import math

def solve():
    """
    This function solves the geometry problem based on the derived value of AE.
    The derivation steps are outlined above. The key results are:
    AC^2 = 120
    cos(angle BAC) = 1/4
    AE = AC * sin(angle ABC) / sin(angle BAC)
    """

    # From detailed analysis, the side length AC=AB=S is found to be S=sqrt(120).
    S_squared = 120
    S = math.sqrt(S_squared)
    
    # cos(alpha) where alpha is angle BAC is 1/4
    cos_alpha = 1/4
    sin_alpha = math.sqrt(1 - cos_alpha**2) # sqrt(15)/4
    
    # In triangle ABC, angle ABC = beta. alpha + 2*beta = 180 degrees.
    # beta = 90 - alpha/2.
    # sin(beta) = sin(90-alpha/2) = cos(alpha/2)
    # cos(alpha/2) = sqrt((1+cos_alpha)/2)
    cos_alpha_over_2 = math.sqrt((1 + cos_alpha) / 2)
    sin_beta = cos_alpha_over_2
    
    # In triangle ABE, it can be shown that angle AEB = alpha.
    # By Law of Sines in triangle ABE: AE/sin(beta) = AB/sin(alpha)
    # AE = AB * sin(beta) / sin(alpha) = S * sin(beta) / sin(alpha)
    AE = S * sin_beta / sin_alpha
    
    # Simplify AE to the form a*sqrt(b)/c
    # AE = sqrt(120) * sqrt(5/8) / (sqrt(15)/4)
    # AE = 2*sqrt(30) * (sqrt(5)/(2*sqrt(2))) / (sqrt(15)/4)
    # AE = sqrt(30) * sqrt(5)/sqrt(2) * 4/sqrt(15)
    # AE = sqrt(15)*sqrt(2) * sqrt(5)/sqrt(2) * 4/sqrt(15)
    # AE = sqrt(15)*sqrt(5) * 4/sqrt(15) = 4*sqrt(5)
    
    a = 4
    b = 5
    c = 1
    
    # We need to find a+b+c
    result = a + b + c
    
    # To be clear in the output, let's print the equation for AE
    print(f"The length of AE is {a} * sqrt({b}) / {c}")
    print(f"a = {a}")
    print(f"b = {b}")
    print(f"c = {c}")
    print(f"The sum a + b + c is {a} + {b} + {c} = {result}")

solve()