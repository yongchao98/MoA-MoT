import math

def solve_magnetization():
    """
    Calculates the minimum magnetization M_z(1) based on the derived formula.
    The minimum occurs at n=1, and the value is -32 / pi^2.
    """
    
    numerator = -32
    denominator_val = math.pi
    denominator_pow = 2
    
    result = numerator / (denominator_val ** denominator_pow)
    
    print(f"The equation for the minimum magnetization is:")
    print(f"M_z(1) = {numerator} / (π^{denominator_pow})")
    print(f"Calculating the value:")
    print(f"M_z(1) = {numerator} / ({denominator_val:.4f}^{denominator_pow}) = {result:.4f}")

solve_magnetization()
