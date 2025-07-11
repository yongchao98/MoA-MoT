import math

# Step 1: Define the given system parameters.
# c: Consistency level of the decentralized identifier resolution mechanism.
c = 0.95
# b: Branching factor of semantic version control.
b = 3

print("This script calculates the theoretical maximum FAIR compliance score (R) for a federated knowledge graph.")
print(f"The calculation is based on a consistency level (c) of {c} and a branching factor (b) of {b}.\n")

# Step 2: Model the maximum theoretical values for each FAIR component.

# Findability (f) is limited by the consistency of identifier resolution.
f = c

# Accessibility (a) is assumed to be perfect (1.0) under "best implementation".
a = 1.0

# Interoperability (i) is inversely proportional to the semantic ambiguity
# introduced by the branching factor 'b'.
i = 1 / b

# Reusability (r) is also hindered by the same semantic ambiguity.
r = 1 / b

# Step 3: Calculate the overall FAIR score (R) as the average of the components.
R = (f + a + i + r) / 4

# Step 4: Display the breakdown and the final result.
print("The FAIR score (R) is the average of its four components: R = (f + a + i + r) / 4")
print("-" * 70)
print(f"Findability (f)      = c               = {f:.4f}")
print(f"Accessibility (a)    = ideal           = {a:.4f}")
print(f"Interoperability (i) = 1/b             = {i:.4f}")
print(f"Reusability (r)      = 1/b             = {r:.4f}")
print("-" * 70)

# Display the final equation with the numbers plugged in.
print("\nFinal Equation:")
print(f"R = ({f} + {a} + {i:.4f} + {r:.4f}) / 4")
print(f"R = {R:.4f}")
