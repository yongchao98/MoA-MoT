import math

# Based on the problem's physical and geometrical setup, the moment the chain
# loses contact with the ground is determined by solving a quadratic equation for s = sin(t),
# where t is the time in seconds.
# The equation is of the form: a*s^2 + b*s + c = 0.

# Define the coefficients of the quadratic equation.
# These coefficients are derived from the parameters given in the problem:
# robot height (h=1m), arm length (r=0.25m), chain length (lc=10m),
# path radius (R=10m), path tilt angle (alpha=30 deg), robot and arm angular velocities (1 rad/s).
a = 1.0
b = -(40 + math.sqrt(3))
c = 39 - 4 * math.sqrt(3)

# Print the equation being solved.
print("To find the time 't', we solve the following quadratic equation for s = sin(t):")
print(f"({a:.2f})s^2 + ({b:.2f})s + ({c:.2f}) = 0\n")

# Calculate the discriminant to find the roots.
discriminant = b**2 - 4 * a * c

# The quadratic formula gives two solutions for s. We need the physically possible one (|s|<=1).
# We select the smaller root using the '-' sign in the numerator, as the other root is > 1.
s = (-b - math.sqrt(discriminant)) / (2 * a)

# The first time t>0 corresponds to the principal value of the arcsin function.
t = math.asin(s)

# Print the results of each calculation step.
print("The coefficients of the equation are:")
print(f"a = {a}")
print(f"b = {b}")
print(f"c = {c}\n")

print(f"Solving for s = sin(t), the valid solution is:")
print(f"sin(t) = {s}\n")

print("The first time 't' when the chain loses contact is:")
print(f"t = arcsin({s})")
print(f"t = {t} seconds")

# Final Answer
# print(f'<<<{{ {t} }}>>>')