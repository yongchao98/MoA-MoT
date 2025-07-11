import math
from fractions import Fraction

# Step 1: Define initial parameters and calculate the area of the first circle.
w0 = 6
h0 = 8
r0 = w0 / 3
area0 = math.pi * r0**2
area0_coeff = Fraction(int(r0**2))

print(f"Initial rectangle width: {w0}, height: {h0}")
print(f"Radius of the first circle (r0 = w0 / 3): {r0}")
print(f"Area of the first circle (A0 = pi * r0^2): {area0_coeff} * pi")
print("-" * 30)

# Step 2: Determine the scaling factor for subsequent rectangles.
# Diagonal length of the initial rectangle.
d0 = math.sqrt(w0**2 + h0**2)
# The ratio of the circle's radius to the half-length of the diagonal.
k = r0 / (d0 / 2)
# The linear scaling factor for the dimensions of the new rectangles.
s = (1 - k) / 2
s_frac = Fraction(s).limit_denominator()

print(f"Length of the diagonal: {d0}")
print(f"Ratio of radius to half-diagonal (k): {k:.2f}")
print(f"Linear scaling factor for new rectangles (s = (1-k)/2): {s_frac}")
print("-" * 30)

# Step 3: Formulate and calculate the common ratio of the geometric series.
# At each step, 4 new circles are added, and their area is scaled by s^2.
common_ratio = 4 * s**2
common_ratio_frac = Fraction(common_ratio).limit_denominator()

print("The total area is the sum of an infinite geometric series.")
print(f"The first term is a = {area0_coeff} * pi.")
print(f"The common ratio is R = 4 * s^2 = 4 * ({s_frac})^2 = {common_ratio_frac}.")
print("-" * 30)

# Step 4: Calculate the sum of the series.
# Sum = a / (1 - R)
total_area_coeff = area0_coeff / (1 - common_ratio_frac)
total_area = float(total_area_coeff) * math.pi

print("The sum of the series is S = a / (1 - R).")
print(f"S = ({area0_coeff} * pi) / (1 - {common_ratio_frac})")
print(f"S = ({area0_coeff} * pi) / ({1 - common_ratio_frac})")
print(f"S = {total_area_coeff} * pi")
print("-" * 30)

# Final Answer
print(f"The final value of the total area is: {total_area}")
