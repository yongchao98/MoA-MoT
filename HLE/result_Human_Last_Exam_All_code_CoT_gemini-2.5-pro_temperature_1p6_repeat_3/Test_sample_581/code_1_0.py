# The best known lower bound for the size of a cap set in dimension 8, r_3(8),
# is established through the product construction, combining the results for
# dimension 7 and dimension 1.

# The size of the largest cap set in dimension 1 is exactly 2.
r3_1_lower_bound = 2

# The best known lower bound for the size of a cap set in dimension 7,
# found by Ferret and Storme (2002), is 248.
r3_7_lower_bound = 248

# The product construction r_3(n+m) >= r_3(n) * r_3(m) gives us the lower bound for n=8.
# We set n=7 and m=1.
r3_8_lower_bound = r3_7_lower_bound * r3_1_lower_bound

# Print the final equation and the result.
print(f"The best known lower bound for the size of a cap set in dimension 8 is found using the product construction with dimensions 7 and 1.")
print(f"The equation for this construction is: {r3_7_lower_bound} * {r3_1_lower_bound} = {r3_8_lower_bound}")
print(f"The resulting lower bound is {r3_8_lower_bound}.")
