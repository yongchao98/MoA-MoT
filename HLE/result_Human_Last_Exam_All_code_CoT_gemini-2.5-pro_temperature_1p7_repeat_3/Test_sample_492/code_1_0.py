import math

# My derivation shows that the average degree of a vertex at time t, lambda(t),
# is given by the formula: lambda(t) = t^2 / 3.
# The giant component emerges when the average degree becomes 1.
# We solve for the critical time c by setting lambda(c) = 1.

# This gives the equation:
# c^2 / 3 = 1

print("The critical condition for the emergence of the giant component is given by the equation:")
print("c^2 / A = B")

denominator = 3
right_hand_side = 1

print(f"In this equation, the value of A is {denominator}.")
print(f"The value of B is {right_hand_side}.")

# To solve for c, we first rearrange the equation to find c^2.
# c^2 = 3 * 1 = 3
c_squared = denominator * right_hand_side
print(f"Multiplying both sides by A gives the equation: c^2 = {c_squared}")

# The final step is to take the square root. Since time cannot be negative,
# we take the positive root.
c = math.sqrt(c_squared)
print(f"The exact value of c is the square root of {c_squared}.")
print(f"The numerical value of c is approximately {c}.")