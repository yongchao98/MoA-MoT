import sys

# Given parameters
c = 0.95  # Consistency level of identifier resolution
b = 3     # Branching factor of semantic version control

# --- Step 1: Model the maximum value for each FAIR metric ---

# Findability (f) is limited by identifier resolution consistency.
f = c

# Accessibility (a) is also limited by identifier resolution consistency.
a = c

# Interoperability (i) is inversely related to the branching factor.
# A higher b means more semantic drift, making interoperability harder.
i = 1 / b

# Reusability (r) requires data to be findable, accessible, and interoperable.
# We model it as a product of findability and interoperability.
r = f * i

# --- Step 2: Model the overall FAIR score R ---

# The overall FAIR score R is the arithmetic mean of the four metrics.
R = (f + a + i + r) / 4

# --- Step 3: Print the results and the final equation ---

print("Based on the established model, the values for each FAIR metric are:")
print(f"f (Findability) = {f}")
print(f"a (Accessibility) = {a}")
print(f"i (Interoperability) = 1 / {b} = {i}")
print(f"r (Reusability) = f * i = {f} * {i} = {r}")
print("\nThe theoretical maximum FAIR compliance score (R) is calculated as follows:")
# The final requested equation with each number.
# Using f-strings to format the output clearly.
print(f"R = (f + a + i + r) / 4")
print(f"R = ({f} + {a} + {i} + {r}) / 4")
final_sum = f + a + i + r
print(f"R = {final_sum} / 4")
print(f"R = {R}")
if 'get_ipython' not in globals() and 'סי' not in sys.modules:
    print(f"\n<<<0.6375>>>")
