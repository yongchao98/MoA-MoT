import math

# Given parameters
# c: consistency level of the decentralized identifier resolution
c = 0.95
# b: branching factor of semantic version control
b = 3

# Step 1: Model Findability (f)
# Findability is limited by the consistency of the identifier resolution.
f = c

# Step 2: Model Accessibility (a)
# Assuming best implementation, accessibility via open protocols is maximal.
a = 1.0

# Step 3: Model Interoperability (i)
# Interoperability is proportional to identifier consistency (c)
# and inversely proportional to the semantic fragmentation (b).
i = c / b

# Step 4: Model Reusability (r)
# Reusability is proportional to identifier consistency (c) for linking/provenance
# and inversely proportional to the complexity of the versioning history (b).
r = c / b

# Step 5: Calculate the final FAIR score R as the arithmetic mean.
R = (f + a + i + r) / 4

# Output the equation with the calculated values
print("The FAIR compliance score (R) is calculated as the mean of its four components.")
print("R = (f + a + i + r) / 4")
print(f"Substituting the modeled values based on c = {c} and b = {b}:")
# Note: Using formatted strings to show the equation clearly.
print(f"R = ({f} + {a} + {i:.4f} + {r:.4f}) / 4")
print(f"R = ({f} + {a} + {c}/{b} + {c}/{b}) / 4")

# Output the final result
print(f"\nThe theoretical maximum value of R is: {R}")
