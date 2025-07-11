import math

# Given values
c1 = 10**4
# The value l1 is given by the expression (1 + 10^5)^5
l1 = (1 + 10**5)**5

# The derived equation for u1 is u1 = (x11 - 1) / c1
# We assume x11 = l1
# So, u1 = (l1 - 1) / c1

# Calculate u1
# Python's integers handle large numbers, so l1 will be calculated exactly.
u1 = (l1 - 1) / c1

# Print the final equation with all numbers
print(f"The equation for u1 is:")
print(f"u1 = (l1 - 1) / c1")
print(f"Substituting the values:")
print(f"u1 = (({1 + 10**5})^5 - 1) / {c1}")
print(f"u1 = ({100001}^5 - 1) / {c1}")
print(f"u1 = ({l1} - 1) / {c1}")
print(f"u1 = {l1 - 1} / {c1}")
print(f"The final value is:")
print(f"u1 = {u1}")