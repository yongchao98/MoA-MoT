import math
from fractions import Fraction

def solve_circle_area_limit():
    """
    Calculates the limit of the total area of circles in the described fractal process.
    """
    # Step 1: Calculate the area of the first circle (the first term, a).
    print("Step 1: Calculate the area of the first circle.")
    
    # Initial rectangle dimensions
    w0 = 6
    
    # The radius of the first circle is 1/3 of the width.
    # We use Fractions for precision.
    r0 = Fraction(w0, 3)
    
    # The area of the first circle is pi * r0^2. We'll track the coefficient of pi.
    # This is the first term 'a' of our geometric series.
    a_coeff = r0**2
    
    print(f"The initial width is w_0 = {w0}.")
    print(f"The radius of the first circle is r_0 = w_0 / 3 = {w0}/3 = {r0}.")
    print(f"The area of the first circle is a = pi * ({r0})^2 = {a_coeff}*pi.")
    print("-" * 50)

    # Step 2: Determine the common ratio 'r' of the geometric series.
    print("Step 2: Determine the common ratio of the series.")
    
    # At each step, new rectangles are formed in the corners. We need to find their size.
    # For a parent rectangle of width w, the new smaller rectangle has width w_new = (3/10)*w.
    dim_scale_factor = Fraction(3, 10)
    
    # The area of a circle scales by the square of the dimension scaling factor.
    area_scale_factor = dim_scale_factor**2
    
    # The number of new circles multiplies by 4 at each step.
    num_multiplier = 4
    
    # The common ratio 'r' is the product of the number multiplier and the area scaling factor.
    common_ratio = num_multiplier * area_scale_factor
    
    print(f"The number of circles added at each step is multiplied by {num_multiplier}.")
    print(f"The area of each individual new circle is scaled by a factor of ({dim_scale_factor})^2 = {area_scale_factor}.")
    print(f"Thus, the total area added at each step is scaled by a common ratio r = {num_multiplier} * {area_scale_factor} = {common_ratio}.")
    print("-" * 50)

    # Step 3: Calculate the sum of the infinite geometric series.
    print("Step 3: Calculate the sum of the infinite geometric series.")
    
    # The total area S is the sum of the series S = a / (1 - r).
    # We calculate the coefficient of pi for the total area.
    total_area_coeff = a_coeff / (1 - common_ratio)
    
    print("The total area S is the sum of the series a + a*r + a*r^2 + ...")
    print("The formula for the sum is S = a / (1 - r).")
    
    # Display the final equation with all numbers.
    a_num = a_coeff.numerator
    a_den = a_coeff.denominator
    r_num = common_ratio.numerator
    r_den = common_ratio.denominator
    one_minus_r = 1 - common_ratio
    
    print("\nFinal Equation:")
    print(f"Total Area = ({a_coeff} * pi) / (1 - {common_ratio})")
    print(f"           = ({a_coeff} * pi) / ({one_minus_r})")
    print(f"           = ({a_num} * {one_minus_r.denominator} * pi) / ({a_den} * {one_minus_r.numerator})")
    print(f"           = ({a_num * one_minus_r.denominator} * pi) / {a_den * one_minus_r.numerator}")
    print(f"           = {total_area_coeff.numerator}*pi / {total_area_coeff.denominator}")

    # Calculate the final numerical value.
    final_value = float(total_area_coeff) * math.pi
    print(f"\nThe exact total area is {total_area_coeff}*pi, which is approximately {final_value:.5f}.")

solve_circle_area_limit()