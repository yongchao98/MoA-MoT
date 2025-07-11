import math

# The problem is to find the constant 'b' in the asymptotic formula
# C(n) ~ b * n^(3/2), where C(n) is the expected cover-and-return time
# for a random walk on a random n-vertex tree.

# Based on the theoretical derivation outlined above, the constant b
# is given by the formula sqrt(2 * pi).

# The final equation is b = sqrt(2 * pi)
# Here are the numbers in the equation:
num1 = 2
pi_val = math.pi

# Now, we calculate the value of b.
b = math.sqrt(num1 * pi_val)

print(f"The equation for the constant b is: b = sqrt({num1} * {pi_val:.5f}...)")
print(f"The exact value of b is sqrt(2*pi).")
print(f"The numerical value of b is: {b}")
