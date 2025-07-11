import math

# This script calculates the constant 'c', the time of emergence of the
# giant component in the described random graph model.

# The analysis shows that the average degree of the graph at time 't'
# is given by the function k(t) = t^2 / 3.
# The giant component emerges when the average degree k(t) equals 1.
# We solve the equation c^2 / 3 = 1 to find the critical time 'c'.

# --- Calculation Steps ---

# The equation to solve:
variable_squared = "c^2"
denominator = 3
target_degree = 1

print("The condition for the emergence of the giant component is that the average degree equals 1.")
print(f"We set the average degree at time 'c' to 1, which gives the equation:")
print(f"{variable_squared} / {denominator} = {target_degree}")

# Solve for c^2
c_squared = target_degree * denominator
print("\nStep 1: Solve for c^2.")
print(f"{variable_squared} = {target_degree} * {denominator}")
print(f"{variable_squared} = {c_squared}")

# Solve for c
c = math.sqrt(c_squared)
print("\nStep 2: Solve for c by taking the square root.")
print(f"c = sqrt({c_squared})")

# Final Answer
print("\nThe exact value of c is the square root of 3.")
print(f"The numerical value of c is: {c}")