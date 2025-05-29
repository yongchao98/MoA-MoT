# Given values
sum_of_edges = 28
sum_of_products = 192

# Calculate a^2 + b^2 + c^2
a2_b2_c2 = sum_of_edges**2 - 2 * sum_of_products

# Calculate the radius r
r = (a2_b2_c2**0.5) / 2

print(r)