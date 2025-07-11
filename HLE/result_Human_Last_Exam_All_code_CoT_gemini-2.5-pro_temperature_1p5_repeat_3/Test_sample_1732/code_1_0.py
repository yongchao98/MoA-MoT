import math

# Case 1 from the derivation leads to a quadratic equation for d.
# The equation is 50d^2 - d - 51 = 0.
# We will solve this equation and apply the condition d > 1.

# Coefficients of the quadratic equation a*x^2 + b*x + c = 0
a = 50
b = -1
c = -51

print("The final equation to solve for the common difference 'd' is:")
# Output each number in the final equation as requested.
print(f"{a}d^2 + ({b})d + ({c}) = 0")
print("")

# Using the quadratic formula to solve for d: d = (-b +/- sqrt(b^2 - 4ac)) / (2a)
discriminant = b**2 - 4 * a * c
sol1 = (-b + math.sqrt(discriminant)) / (2 * a)
sol2 = (-b - math.sqrt(discriminant)) / (2 * a)

print(f"The possible values for d are {sol1} and {sol2}.")
print("")

# Check against the given condition d > 1
print("Applying the condition that d > 1:")
final_d = None
if sol1 > 1:
    final_d = sol1
    print(f"The solution {sol1} is greater than 1.")
if sol2 > 1:
    final_d = sol2
    print(f"The solution {sol2} is greater than 1.")

if final_d is not None:
    print(f"\nSince only {final_d} satisfies the condition d > 1, this is our answer.")
    print(f"The value of d is {final_d}.")
else:
    print("No solution from this case satisfies d > 1.")
