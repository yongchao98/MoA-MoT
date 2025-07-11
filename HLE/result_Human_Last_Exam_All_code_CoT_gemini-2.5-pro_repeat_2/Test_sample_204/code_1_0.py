import math

# The Hopf charge is the linking number of the preimages of two points.
# Let's choose the South Pole (0,0,-1) and an equatorial point (1,0,0).

# Preimage of South Pole: nz = cos(G) = -1 => G = pi
# G = pi * exp(-10*r2) => exp(-10*r2) = 1 => r2 = 0
# r2^2 = (x*x+y*y-0.5)^2 + z*z = 0
# This means z=0 and x*x+y*y-0.5 = 0.
# So, C1 is the circle: x^2 + y^2 = 0.5 in the z=0 plane.
radius_C1_sq = 0.5
print(f"The first curve, C1 (preimage of the South Pole), is a circle in the z=0 plane with equation x^2 + y^2 = {radius_C1_sq}.")

# Preimage of (1,0,0): nz = cos(G) = 0 => G = pi/2
# G = pi * exp(-10*r2) => exp(-10*r2) = 0.5 => r2 = -ln(0.5)/10 = ln(2)/10
# Also nx = sin(G)cos(f) = 1. Since G=pi/2, sin(G)=1, so cos(f)=1.
# This means f=0, which corresponds to y=0 and x>0.
# So, C2 is the curve (x^2-0.5)^2 + z^2 = (ln(2)/10)^2 in the y=0, x>0 half-plane.
Rh = math.log(2) / 10
print(f"The second curve, C2 (preimage of (1,0,0)), is a loop in the y=0 plane with equation (x^2 - 0.5)^2 + z^2 = {Rh**2:.6f}.")

# To find the linking number, we check where C2 pierces the disk bounded by C1.
# The disk D1 is x^2+y^2 <= 0.5 in the z=0 plane.
# C2 pierces the z=0 plane when (x^2-0.5)^2 = Rh^2.
# x^2 - 0.5 = +/- Rh
# x1_sq = 0.5 - Rh
# x2_sq = 0.5 + Rh
x1_sq = 0.5 - Rh
x2_sq = 0.5 + Rh

print(f"\nC2 pierces the z=0 plane at two points where x^2 has values:")
print(f"x1^2 = 0.5 - ln(2)/10 = {x1_sq:.6f}")
print(f"x2^2 = 0.5 + ln(2)/10 = {x2_sq:.6f}")

# Check if these piercing points are inside the disk D1 (x^2 <= 0.5).
is_inside_1 = x1_sq <= radius_C1_sq
is_inside_2 = x2_sq <= radius_C1_sq

print(f"\nIs the first piercing point (x1^2={x1_sq:.6f}) inside the disk (x^2<={radius_C1_sq})? {is_inside_1}")
print(f"Is the second piercing point (x2^2={x2_sq:.6f}) inside the disk (x^2<={radius_C1_sq})? {is_inside_2}")

# Since one point is inside and one is outside, the curve C2 pierces the disk D1 exactly once.
hopf_charge = 1
print(f"\nSince the curve C2 pierces the disk bounded by C1 exactly once, the linking number is 1.")
print(f"The Hopf charge is {hopf_charge}.")
