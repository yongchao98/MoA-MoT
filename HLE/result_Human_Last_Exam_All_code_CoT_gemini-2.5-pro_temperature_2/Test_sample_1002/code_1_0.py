# The problem is stated for a general integer k >= 2.
# We will represent the solution symbolically in terms of 'k'.
k_symbol = 'k'

# From the mathematical derivation based on extremal graph theory, the limit
# is a rational expression in terms of k.
numerator_expr = k_symbol
denominator_expr = f"({k_symbol} + 1)"

# We print the components of the final expression as requested.
print("The limit is a symbolic expression depending on k.")
print(f"The numerator of the final fraction is: {numerator_expr}")
print(f"The denominator of the final fraction is: {denominator_expr}")
print(f"Therefore, the final limit is the ratio: {numerator_expr} / {denominator_expr}")
