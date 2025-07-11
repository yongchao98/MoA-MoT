# Given system parameters
c = 0.95  # Consistency level of identifier resolution
b = 3     # Branching factor of semantic version control

# Step 1: Model the individual FAIR metrics based on the parameters.
# Findability is limited by identifier consistency.
f = c
# Accessibility is also limited by identifier consistency.
a = c
# Interoperability is inversely affected by semantic branching.
i = 1 / b
# Reusability requires finding the data and it being interoperable.
r = f * i

# Step 2: Define the overall FAIR score R as the average of the four metrics.
R = (f + a + i + r) / 4

# Step 3: Print the final calculation and the result.
# The f-string formats the numbers to a reasonable number of decimal places for clarity.
print(f"The FAIR compliance score R is calculated as (f + a + i + r) / 4.")
print(f"R = ({f:.2f} + {a:.2f} + {i:.3f} + {r:.3f}) / 4")
print(f"R = {f + a + i + r:.4f} / 4")
print(f"The theoretical maximum value of R is: {R:.4f}")
