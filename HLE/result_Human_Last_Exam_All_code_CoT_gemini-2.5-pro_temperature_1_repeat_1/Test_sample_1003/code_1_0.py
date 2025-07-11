import math

# The problem asks for the value of the ratio R = (1 - cos(θ'_14)) / (1 - cos(θ'_34)).
# As explained in the plan, this complex ratio simplifies to a much more manageable expression
# involving the known angles between stars S1, S2, and S3:
# R = (1 - cos(θ'_12)) / (1 - cos(θ'_13))

# We are given the angles in the second frame:
# θ'_12 = π/2 (right angle)
# θ'_13 = 3π/4

# Let's define the numbers that will be part of the final equation.
a = 1.0
b = math.cos(math.pi / 2)  # cos(π/2) = 0
c = 1.0
d = math.cos(3 * math.pi / 4) # cos(3π/4) = -sqrt(2)/2

# The final equation to solve is of the form (a - b) / (c - d)
print("The final equation for the ratio is: (a - b) / (c - d)")
print("\nThe numbers in this final equation are:")
print(f"a = {a}")
print(f"b = cos(π/2) = {b}")
print(f"c = {c}")
print(f"d = cos(3π/4) = {d}")

# Now, we calculate the final value by substituting these numbers.
result = (a - b) / (c - d)

print(f"\nThus, the calculation is: ( {a} - {b} ) / ( {c} - ({d}) )")
print(f"Result: {result}")