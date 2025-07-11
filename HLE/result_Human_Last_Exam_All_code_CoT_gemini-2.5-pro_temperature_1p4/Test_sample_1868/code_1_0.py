# Define the given parameters
c = 0.95  # Consistency level of identifier resolution
b = 3     # Branching factor of semantic version control

# Model the FAIR components based on the problem description
# Findability and Accessibility are limited by the consistency 'c'.
f = c
a = c

# Interoperability is inversely affected by the branching factor 'b'.
i = 1 / b

# Reusability is assumed to be perfect (1.0) under the "best implementation" scenario.
r = 1.0

# Calculate the final FAIR compliance score R as the arithmetic mean of the components.
R = (f + a + i + r) / 4

# Output the final equation and result, showing each number in the equation.
# We format the numbers to make the equation clear.
print(f"The theoretical maximum FAIR score R is calculated as follows:")
print(f"R = (f + a + i + r) / 4")
print(f"R = ({f:.2f} + {a:.2f} + {i:.2f} + {r:.2f}) / 4")
final_value_str = f"({f} + {a} + {1/b:.3f} + {r})"
print(f"R = ({f} + {a} + 1/{b} + {r}) / 4")
print(f"R = {R}")