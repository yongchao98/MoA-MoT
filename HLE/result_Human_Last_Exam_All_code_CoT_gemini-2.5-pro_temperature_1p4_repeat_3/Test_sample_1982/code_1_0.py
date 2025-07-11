import math

# Step 1: Define the given constants.
# We will use Python's built-in arbitrary-precision integers for all calculations
# to ensure accuracy with very large numbers.
c1 = 10**4
r = 10**5

# Step 2: Define and calculate alpha1.
# The formula is alpha1 = (1 + 10^5)^6 * (1 - 10^5 + 10^10).
# A simplified but equivalent form is alpha1 = (1 + 10^5)^5 * (1 + (10^5)^3),
# which we can derive from the sum of cubes formula: (a+b)(a^2-ab+b^2) = a^3+b^3.
# Let's use the original formula as given in the problem.
alpha1 = (1 + r)**6 * (1 - r + r**2)

# Step 3: Set x11 based on our assumption.
# The formula for u1 depends on x11, which we assume to be alpha1.
x11 = alpha1

# Step 4: From the matrix equation, we derive the formula for u1:
# x11 = 1 + c1 * u1  =>  u1 = (x11 - 1) / c1
# We can confirm that (x11 - 1) is exactly divisible by c1.
u1 = (x11 - 1) // c1

# Step 5: Print the final equation with the computed numerical values.
print("Derived formula for the control u1:")
print("u1 = (x11 - 1) / c1")
print("\nAssuming x11 = alpha1, we substitute the values:")
print(f"c1 = {c1}")
print(f"alpha1 = {alpha1}")
print("\nThe final equation with numbers is:")
print(f"u1 = ({alpha1} - 1) / {c1}")
print("\nResult:")
print(f"u1 = {u1}")