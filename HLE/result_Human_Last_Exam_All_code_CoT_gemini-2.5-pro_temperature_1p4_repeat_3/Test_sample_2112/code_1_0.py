import math

def solve_r0():
    """
    This function calculates the value of r0 based on the derived equation.
    The equation for r0 is: r0 = (226 + 49 * sqrt(2)) / 17
    """
    
    numerator_const = 226
    sqrt_coeff = 49
    sqrt_val = math.sqrt(2)
    denominator = 17
    
    r0 = (numerator_const + sqrt_coeff * sqrt_val) / denominator
    
    print("The radial distance r_0 is determined by the equation:")
    print(f"r_0 = ({numerator_const} + {sqrt_coeff} * sqrt({int(sqrt_val**2)})) / {denominator}")
    print("\nThe numerical value is:")
    print(r0)

solve_r0()
