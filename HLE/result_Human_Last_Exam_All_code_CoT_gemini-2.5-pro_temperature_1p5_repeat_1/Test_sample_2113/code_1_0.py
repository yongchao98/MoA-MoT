import math

# Plan:
# 1. Based on the harmonic potential, the trajectories are simple harmonic motions.
#    Given z(0) = sqrt(2) and z1(0) = 1, and assuming release from rest (z'(0)=0),
#    we get z(t) = sqrt(2)*cos(t) and z1(t) = cos(t).
# 2. We solve the integral equation for y(t) using z(t) = sqrt(2)*cos(t).
#    This yields y(t) = sin(2t) * (cos(2t))^(-3/2).
# 3. We calculate (z1(pi/8) / y(pi/8))^2.

# Value of t
t = math.pi / 8
sqrt2 = math.sqrt(2)

# Calculation of the numerator: z1(pi/8)^2
# z1(pi/8) = cos(pi/8)
# We need cos(pi/8)^2. Using the half-angle formula cos(x)^2 = (1 + cos(2x))/2:
# cos(pi/8)^2 = (1 + cos(pi/4))/2 = (1 + 1/sqrt(2))/2
z1_squared = (1 + math.cos(math.pi / 4)) / 2

# Calculation of the denominator: y(pi/8)^2
# y(pi/8) = sin(2*pi/8) * cos(2*pi/8)^(-3/2)
# y(pi/8) = sin(pi/4) * cos(pi/4)^(-3/2)
# sin(pi/4) = 1/sqrt(2), cos(pi/4) = 1/sqrt(2)
# y(pi/8) = (1/sqrt(2)) * (1/sqrt(2))^(-3/2) = (2^(-0.5)) * (2^(-0.5))^(-1.5) = 2^(-0.5) * 2^(0.75) = 2^0.25
# y(pi/8)^2 = (2^0.25)^2 = 2^0.5 = sqrt(2)
y_squared = math.sqrt(2)

# Final calculation
final_result = z1_squared / y_squared

print("The problem simplifies by identifying the particle motion as simple harmonic.\n")
print(f"We need to compute (z1(\u03C0/8) / y(\u03C0/8))^2.\n")

print("Step 1: Calculate the numerator, z1(\u03C0/8)^2")
print("z1(\u03C0/8)^2 = cos(\u03C0/8)^2")
print(f"Using the half-angle identity, cos(x)^2 = (1 + cos(2x))/2:")
print(f"cos(\u03C0/8)^2 = (1 + cos(\u03C0/4))/2 = (1 + 1/\u221A2)/2")
print(f"Symbolically, z1(\u03C0/8)^2 = (\u221A2 + 1)/(2*\u221A2)")
print(f"The numerical value is: {z1_squared:.7f}\n")


print("Step 2: Calculate the denominator, y(\u03C0/8)^2")
print("From the integral equation, y(t) = sin(2t) * (cos(2t))^(-3/2).")
print("y(\u03C0/8) = sin(\u03C0/4) * (cos(\u03C0/4))^(-3/2) = (1/\u221A2) * (1/\u221A2)^(-3/2) = 2^(1/4)")
print(f"Therefore, y(\u03C0/8)^2 = (2^(1/4))^2 = \u221A2")
print(f"The numerical value is: {y_squared:.7f}\n")


print("Step 3: Compute the final ratio")
print("(z1(\u03C0/8) / y(\u03C0/8))^2 = z1(\u03C0/8)^2 / y(\u03C0/8)^2")
# The final equation printed with symbolic values
print("= [(\u221A2 + 1)/(2*\u221A2)] / [\u221A2]")
print("= (\u221A2 + 1) / (2 * \u221A2 * \u221A2)")
print("= (\u221A2 + 1) / 4\n")

print("The final numerical result is:")
print(f"({sqrt2:.7f} + 1) / 4 = {final_result:.7f}")
