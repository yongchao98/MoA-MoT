import math
from fractions import Fraction

def solve_circle_area():
    """
    Calculates the limit of the total area of all circles in the fractal R_n.
    """
    # Step 1: Define initial parameters and calculate the area of the first circle (A1).
    w0 = 6
    r1 = Fraction(w0, 3)
    # Area of the first circle, A1 = pi * r1^2. We handle pi symbolically.
    a_numeric = r1**2
    
    print(f"The initial rectangle has width w_0 = {w0}.")
    print(f"The radius of the first circle is r_1 = w_0 / 3 = {w0} / 3 = {r1}.")
    print(f"The area of the first circle is A_1 = pi * ({r1})^2 = {a_numeric}*pi.")
    print("-" * 30)

    # Step 2: Determine the common ratio 'q' of the geometric series of areas.
    # The linear dimension scaling factor is 3/10.
    scaling_factor_dim = Fraction(3, 10)
    # Area scales by the square of the linear scaling factor.
    scaling_factor_area = scaling_factor_dim**2
    # At each step, 4 new circles are added for each previous circle.
    num_new_circles = 4
    q = num_new_circles * scaling_factor_area

    print("At each step, 4 new smaller rectangles are created.")
    print(f"The linear scaling factor for the new rectangles is {scaling_factor_dim}.")
    print(f"The area scaling factor for the new circles is ({scaling_factor_dim})^2 = {scaling_factor_area}.")
    print(f"The common ratio for the total area added at each step is q = (num new circles) * (area scaling factor) = {num_new_circles} * {scaling_factor_area} = {q}.")
    print("-" * 30)

    # Step 3: Sum the infinite geometric series S = a / (1 - q).
    # 'a' is the first term, which is the area of the first circle.
    a = a_numeric
    denominator = 1 - q

    # Calculate final result
    total_area_numeric = a / denominator
    
    # Present the final calculation step-by-step
    print("The total area is the sum of an infinite geometric series S = a / (1 - q).")
    print(f"Here, a = A_1 = {a}*pi and q = {q}.")
    print("The final equation is:")
    # Output each number in the final equation
    print(f"S = ({a} * pi) / (1 - {q.numerator}/{q.denominator})")
    print(f"S = ({a} * pi) / ({denominator.numerator}/{denominator.denominator})")
    print(f"S = ({a} * pi) * ({denominator.denominator}/{denominator.numerator})")
    final_num = a.numerator * denominator.denominator
    final_den = a.denominator * denominator.numerator
    result_frac = Fraction(final_num, final_den)
    print(f"S = ({final_num}*pi) / {final_den} = {result_frac.numerator}/{result_frac.denominator}*pi")

solve_circle_area()
print("<<<25/4*pi>>>")