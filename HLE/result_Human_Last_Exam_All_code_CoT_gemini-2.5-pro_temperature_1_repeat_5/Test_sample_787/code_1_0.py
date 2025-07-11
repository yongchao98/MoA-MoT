# The reasoning above shows that the limit is the integer L = 2^10 * 3^2 * 5^1.
# This code calculates this value.

# Exponents
a = 10
b = 2
c = 1

# Base values
base_2 = 2
base_3 = 3
base_5 = 5

# Calculation
val_2 = base_2**a
val_3 = base_3**b
val_5 = base_5**c

limit = val_2 * val_3 * val_5

print(f"{base_2}^{a} * {base_3}^{b} * {base_5}^{c} = {val_2} * {val_3} * {val_5} = {limit}")