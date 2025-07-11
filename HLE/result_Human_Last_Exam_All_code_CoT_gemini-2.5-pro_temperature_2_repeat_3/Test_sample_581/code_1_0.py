# The dimension of the vector space over F_3
n = 8

# The best known lower bound for the size of a cap set in dimension 8.
# This value is a result of extensive mathematical research and computer searches,
# rather than a simple formula.
best_known_lower_bound = 496

# The final equation shows the relationship between the maximum cap set size
# in dimension 8, cap(8), and its known lower bound.
# Here we print each number in the final equation.
print(f"cap({n}) >= {best_known_lower_bound}")