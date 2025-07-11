import math
from fractions import Fraction

# Step 1: Define initial parameters and calculate the area of the first circle (the first term 'a').
w0 = 6
h0 = 8
r1 = w0 / 3
first_term_val = math.pi * r1**2
print(f"The first term of the series is the area of the initial circle.")
print(f"a = pi * (width / 3)^2 = pi * ({w0}/3)^2 = {int(first_term_val/math.pi)}*pi")
print("-" * 30)

# Step 2: Determine the common ratio 'r' of the geometric series.
# The number of circles multiplies by 4 at each step.
# The area of each new circle is scaled by a factor k_A.
# The width scaling factor k_w = 3/10.
# The area scaling factor k_A = k_w^2 = (3/10)^2 = 9/100.
k_A = (3/10)**2
# The common ratio r = (number of new circles factor) * (area scaling factor).
common_ratio_val = 4 * k_A
r_frac = Fraction(common_ratio_val).limit_denominator()
print("The common ratio 'r' is the product of the multiplication factor for circles (4) and the area scaling factor.")
print(f"r = 4 * (3/10)^2 = 4 * (9/100) = 36/100 = {r_frac.numerator}/{r_frac.denominator}")
print("-" * 30)

# Step 3: Calculate the sum of the infinite geometric series S = a / (1 - r).
total_area = first_term_val / (1 - common_ratio_val)
a_frac_pi = Fraction(first_term_val/math.pi).limit_denominator()
denominator_frac = 1 - r_frac
sum_frac_pi = a_frac_pi / denominator_frac

print("The total area is the sum of the infinite geometric series S = a / (1 - r).")
print("Here is the calculation with the determined values:")
print(f"S = ({a_frac_pi.numerator}*pi) / (1 - {r_frac.numerator}/{r_frac.denominator})")
print(f"S = ({a_frac_pi.numerator}*pi) / ({denominator_frac.numerator}/{denominator_frac.denominator})")
print(f"S = ({a_frac_pi.numerator * denominator_frac.denominator}*pi) / {denominator_frac.numerator}")
print(f"S = {sum_frac_pi.numerator}*pi / {sum_frac_pi.denominator}")
print("-" * 30)

print(f"The final value of the limit is {sum_frac_pi.numerator}*pi/{sum_frac_pi.denominator}, which is approximately {total_area:.4f}.")
