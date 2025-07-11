import math
from fractions import Fraction

# Step 1: Define initial parameters and calculate the area of the first circle.
w0 = 6
h0 = 8

# The radius of the first circle is 1/3 of the width.
r0 = Fraction(w0, 3)
# The area of the first circle is pi * r0^2. This is the first term 'a' of our series.
# We will calculate the coefficient of pi.
a_coeff = r0**2

# Step 2: Determine the common ratio 'r' of the series.
# The ratio of the width of a new rectangle to the parent rectangle.
# As derived in the plan, w_new = (3/10) * w_old.
width_ratio = Fraction(3, 10)

# The ratio of the area of a new circle to the parent circle's area is width_ratio^2.
area_ratio = width_ratio**2

# At each step, 4 new circles are created for each existing one.
num_circles_multiplier = 4

# The common ratio 'r' is the product of the multiplier and the area ratio.
r = num_circles_multiplier * area_ratio

# Step 3: Calculate the sum of the infinite geometric series.
# The formula is S = a / (1 - r).
# We calculate the coefficient of pi for the final sum.
sum_coeff = a_coeff / (1 - r)

# Step 4: Print the final equation and the result.
print("This problem can be solved by finding the sum of an infinite geometric series S = a / (1 - r).")
print("\n1. Find the first term 'a' (area of the first circle):")
print(f"   Initial width w₀ = {w0}")
print(f"   Initial radius r₀ = w₀ / 3 = {w0}/3 = {r0.numerator}/{r0.denominator}")
print(f"   Area 'a' = π * r₀² = π * ({r0.numerator}/{r0.denominator})² = {a_coeff.numerator}/{a_coeff.denominator} * π")

print("\n2. Find the common ratio 'r':")
print(f"   At each step, we add 4 new circles for each prior one.")
print(f"   The scaling factor for the width of new rectangles is {width_ratio.numerator}/{width_ratio.denominator}.")
print(f"   The scaling factor for the area of a single new circle is ({width_ratio.numerator}/{width_ratio.denominator})² = {area_ratio.numerator}/{area_ratio.denominator}.")
print(f"   The common ratio 'r' = (Number of new circles) * (Area scaling factor) = {num_circles_multiplier} * {area_ratio.numerator}/{area_ratio.denominator} = {r.numerator}/{r.denominator}.")

print("\n3. Calculate the total sum S:")
print(f"   S = a / (1 - r)")
print(f"   S = ({a_coeff.numerator}/{a_coeff.denominator} * π) / (1 - {r.numerator}/{r.denominator})")
print(f"   S = ({a_coeff.numerator}/{a_coeff.denominator} * π) / ({(1-r).numerator}/{(1-r).denominator})")
print(f"   S = ({a_coeff.numerator}/{a_coeff.denominator}) * ({(1-r).denominator}/{(1-r).numerator}) * π")
print(f"   S = {sum_coeff.numerator}/{sum_coeff.denominator} * π")

# Final numerical answer
final_value = float(sum_coeff) * math.pi
print(f"\nThe numerical value is approximately: {final_value}")

# The format requested for the final answer
# Using f-string to format the output string
final_answer_symbolic = f"{sum_coeff.numerator}*pi/{sum_coeff.denominator}"
# The problem asks for the value.
# So we output the numerical value. 25/4 * pi = 6.25 * pi approx 19.635
# Final answer as a symbolic fraction of pi is 25*pi/4
# Final answer as a decimal is approx 19.634954084936208
# The format requested seems to be decimal
print(f"\nFinal Answer: {final_value}")