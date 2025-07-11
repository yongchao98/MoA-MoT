# 1. Define the given constants
c1 = 10**4

# 2. Calculate the value of x11
# The problem provides l1 and alpha1, from which we deduce x11 = alpha1 / l1.
# l1 = (1 + 10**5)**5
# alpha1 = (1 + 10**5)**6 * (1 - 10**5 + 10**10)
# x11 = alpha1 / l1 leads to:
# x11 = (1 + 10**5) * (1 - 1*10**5 + (10**5)**2)
# This simplifies using the sum of cubes identity: a^3 + b^3 = (a+b)(a^2-ab+b^2)
# With a = 1 and b = 10**5, we have:
x11 = 1**3 + (10**5)**3

# 3. Calculate the control u1
# From the matrix equation, we derived the relation: x11 = 1 + c1 * u1
# Solving for u1 gives: u1 = (x11 - 1) / c1
u1 = (x11 - 1) // c1

# 4. Print the final verification equation with each number displayed
print(f"The derived equation relating x11, c1, and u1 is: x11 = 1 + c1 * u1")
print("Substituting the numerical values, we get:")
print(f"{int(x11)} = 1 + {int(c1)} * {int(u1)}")

# 5. Print the final result for u1
print(f"\nTherefore, the value of the control u1 is:")
print(int(u1))