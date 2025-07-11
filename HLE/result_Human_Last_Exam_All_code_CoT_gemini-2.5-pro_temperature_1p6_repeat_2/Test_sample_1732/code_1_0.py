import math

# Based on the derivation, the problem reduces to solving a quadratic equation for d.
# The case a1 = 2d yields solutions d=1 and d=-50/51, neither of which satisfy d > 1.
# The case a1 = d yields the quadratic equation 50d^2 - d - 51 = 0.
# We solve this equation for d.

# Coefficients of the quadratic equation ax^2 + bx + c = 0
a = 50
b = -1
c = -51

print(f"The final equation to solve for d is: {a} * d**2 + {b} * d + {c} = 0")

# Calculate the discriminant
discriminant = b**2 - 4 * a * c

# Check if the discriminant is non-negative for real solutions
if discriminant >= 0:
    # Calculate the two possible solutions for d
    d1 = (-b + math.sqrt(discriminant)) / (2 * a)
    d2 = (-b - math.sqrt(discriminant)) / (2 * a)

    print(f"The potential solutions for d are {d1} and {d2}.")

    # The problem states that d > 1. We find the solution that meets this condition.
    print("\nChecking the condition d > 1...")
    final_d = None
    if d1 > 1:
        final_d = d1
        print(f"The value d = {d1} satisfies the condition d > 1.")
    elif d2 > 1:
        final_d = d2
        print(f"The value d = {d2} satisfies the condition d > 1.")
    else:
         print("Neither solution satisfies the condition d > 1 for this case.")
         
    if final_d is not None:
        print(f"\nThe value of d is {final_d}")

else:
    print("The equation has no real solutions for d.")
