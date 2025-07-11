import math

def get_area_formula():
    """
    Calculates the coefficients of the area formula and prints it.
    The area of the triangle T(t) is A(t) = C1 + C2 * t^2.
    """
    # Given constants from the problem
    R = 10  # Radius of the circumscribed circle
    v = 1   # Speed of the vertices along the hexagon's sides

    # The square of the distance from the center to a vertex is d_v^2 = (3*R^2/4) + t^2
    # The constant part of this expression is:
    constant_term_in_distance_sq = (3 * R**2) / 4

    # The area formula is A(t) = (3 * sqrt(3) / 4) * d_v^2
    # A(t) = (3 * sqrt(3) / 4) * ( (3 * R^2 / 4) + t^2 )
    # Let's substitute R=10 to find the numerical constant.
    # A(t) = (3 * sqrt(3) / 4) * (75 + t^2)

    # Let's define the numbers for the final equation strings
    coeff_num = 3
    coeff_den = 4
    const_term = int(constant_term_in_distance_sq)
    
    # Use the Unicode symbol for the square root for better readability
    sqrt_symbol = "\u221A"

    print("The area of the triangle T(t) as a function of time t is given by the formula:")
    print(f"Area(t) = ({coeff_num}{sqrt_symbol}3 / {coeff_den}) * ({const_term} + t^2)")
    
    # Expanded form: Area(t) = (225*sqrt(3)/4) + (3*sqrt(3)/4)*t^2
    expanded_const_num = coeff_num * const_term
    
    print("\nIn expanded form, this is:")
    print(f"Area(t) = ({expanded_const_num}{sqrt_symbol}3 / {coeff_den}) + ({coeff_num}{sqrt_symbol}3 / {coeff_den}) * t^2")

if __name__ == '__main__':
    get_area_formula()