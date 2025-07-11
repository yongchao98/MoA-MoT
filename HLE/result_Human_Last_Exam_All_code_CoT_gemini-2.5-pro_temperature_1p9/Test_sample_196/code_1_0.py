import math

def solve_volume():
    """
    Calculates the volume of the space enclosed by the cone S1 and the ellipsoid S2.
    The method is based on the finding that the ellipsoid is entirely contained within the cone.
    """

    # The equation for the ellipsoid S2 is x^2/3 + y^2/4 + z^2/3 = 1.
    # From this, we identify the squared values of the semi-axes.
    a_sq = 3
    b_sq = 4
    c_sq = 3

    # The semi-axes are a, b, c.
    a = math.sqrt(a_sq)
    b = math.sqrt(b_sq)
    c = math.sqrt(c_sq)

    # The volume of an ellipsoid is given by the formula V = (4/3) * pi * a * b * c.
    # The numbers to be used in the final equation output are:
    num_4 = 4
    num_3 = 3
    pi_symbol_str = "pi"
    
    # Calculate the volume
    volume = (num_4/num_3) * math.pi * a * b * c
    
    print("The volume is calculated using the ellipsoid volume formula: V = (4/3) * pi * a * b * c")
    print(f"The semi-axes are a = sqrt({a_sq}), b = {int(b)}, c = sqrt({c_sq}).")
    
    print("\nThe final equation for the volume with the values substituted is:")
    # Print the equation step-by-step
    # Step 1: Show the full formula with numbers
    print(f"V = ({num_4}/{num_3}) * {pi_symbol_str} * sqrt({a_sq}) * {int(b)} * sqrt({c_sq})")
    
    # Step 2: Simplify the product of the semi-axes
    product_of_axes_int = int(a * c)
    print(f"V = ({num_4}/{num_3}) * {pi_symbol_str} * {product_of_axes_int} * {int(b)}")
    
    # Step 3: Simplify the numerator of the fraction
    numerator = num_4 * product_of_axes_int * int(b)
    print(f"V = ({numerator}/{num_3}) * {pi_symbol_str}")
    
    # Step 4: Final symbolic result
    final_coefficient = numerator // num_3
    print(f"V = {final_coefficient} * {pi_symbol_str}")
    
    print(f"\nThe final calculated volume is {volume}")

solve_volume()