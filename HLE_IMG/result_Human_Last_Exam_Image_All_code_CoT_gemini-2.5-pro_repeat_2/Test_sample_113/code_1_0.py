import math
from fractions import Fraction

def solve_geometric_series():
    """
    Calculates the limit of the total area of all circles in the described fractal process.
    """
    # Step 1: Define initial parameters for the first rectangle R1.
    w0 = 6
    h0 = 8

    print("Step 1: Analyze the initial figure R1")
    print(f"Initial rectangle width, w0 = {w0}")
    print(f"Initial rectangle height, h0 = {h0}")

    # Calculate the radius and area of the first circle.
    # The radius r1 is 1/3 of the width w0.
    r1 = Fraction(w0, 3)
    # The area A1 is pi * r1^2. This is the first term 'a' of our geometric series.
    # We will calculate the coefficient of pi.
    A1_coeff = r1**2
    print(f"Radius of the first circle, r1 = w0 / 3 = {w0}/3 = {r1}")
    print(f"Area of the first circle, A1 = pi * r1^2 = pi * {r1}^2 = {A1_coeff} * pi")
    print("-" * 50)

    # Step 2: Calculate the common ratio of the geometric series.
    print("Step 2: Determine the common ratio of the series")
    # First, find the scaling factor for the rectangles' dimensions.
    # We need the intersection point of the diagonal and the first circle.
    diag_length = math.sqrt(w0**2 + h0**2)
    
    # x-coordinate of the intersection point in the first quadrant
    x_int = Fraction(w0, diag_length) * r1
    
    # The new smaller rectangle's width w1
    w1 = Fraction(w0, 2) - x_int
    
    # The linear scaling factor 's' is the ratio of the new width to the old width.
    s = w1 / w0
    
    print(f"The diagonal's equation is y = (h0/w0)x and the circle's is x^2 + y^2 = r1^2.")
    print(f"Solving for the intersection gives a new rectangle width w1 = {w1.limit_denominator()}.")
    print(f"The linear scaling factor for dimensions, s = w1 / w0 = {s.limit_denominator()}.")

    # The common ratio 'r' of the geometric series of total areas is 4 * s^2.
    # This is because at each step, we create 4 new rectangles, and their areas are scaled by s^2.
    common_ratio = 4 * s**2
    print(f"The common ratio of the area series, r = 4 * s^2 = 4 * ({s.limit_denominator()})^2 = {common_ratio.limit_denominator()}.")
    print("-" * 50)
    
    # Step 3: Calculate the sum of the infinite geometric series.
    print("Step 3: Sum the infinite geometric series")
    # The sum is S = a / (1 - r), where a = A1_coeff * pi.
    # We calculate the coefficient part first.
    sum_coeff = A1_coeff / (1 - common_ratio)
    
    print("The total area S is the sum of the series a + ar + ar^2 + ... which is S = a / (1 - r).")
    print(f"The first term a = {A1_coeff}*pi and the ratio r = {common_ratio.limit_denominator()}.")
    
    print("\nFinal Equation:")
    print(f"S_infinity = ({A1_coeff} * pi) / (1 - {common_ratio.limit_denominator()})")
    print(f"S_infinity = ({A1_coeff} * pi) / ({(1 - common_ratio).limit_denominator()})")
    print(f"S_infinity = {sum_coeff.limit_denominator()} * pi")
    
    final_value = float(sum_coeff * math.pi)
    print(f"\nThe final numerical value is approximately {final_value:.4f}")

solve_geometric_series()