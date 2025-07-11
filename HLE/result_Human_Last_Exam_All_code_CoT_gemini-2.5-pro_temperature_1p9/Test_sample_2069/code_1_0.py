# Define symbolic placeholders for the bounds l and u.
l = 'l'
u = 'u'

# The first new inequality is y >= l * (1 - b).
# After distributing the brackets, it becomes y >= l - l*b.
inequality1 = f"y >= {l} - {l}*b"

# The second new inequality is y >= x - (u - l) * b.
# (u - l) is a standard choice for the Big-M value here.
# After distributing the brackets, it becomes y >= x - u*b + l*b.
inequality2 = f"y >= x - {u}*b + {l}*b"

# Print the final two inequalities in the required format.
print(f"{inequality1}, {inequality2}")