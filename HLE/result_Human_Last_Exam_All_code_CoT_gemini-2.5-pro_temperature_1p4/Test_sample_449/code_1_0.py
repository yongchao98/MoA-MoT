import math

# Define the starting point coordinates
x = 3000
y = 4000

# Calculate the squared distance from the origin
r_squared = x**2 + y**2

# The probability is given by the formula 1 / (pi^2 * r^2)
# We will show the components of the formula before the final calculation.
pi_squared = math.pi**2
probability = 1 / (pi_squared * r_squared)

print("The problem is to find the probability that a 2D random walk conditioned to avoid the origin, starting far away, never hits the origin's neighbors.")
print("This probability can be approximated by the formula: P = 1 / (pi^2 * |x0|^2)")
print(f"The starting point is x0 = ({x}, {y})")
print(f"The squared distance from the origin is |x0|^2 = {x}^2 + {y}^2 = {int(r_squared)}")
print(f"The value of pi^2 is approximately {pi_squared:.4f}")
print(f"The final equation is P = 1 / ({pi_squared:.4f} * {int(r_squared)})")
print(f"The calculated probability is {probability}")
print(f"The approximate answer with two significant digits is {probability:.2g}")
