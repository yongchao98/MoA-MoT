# Define the given parameters
# c: consistency level of the decentralized identifier resolution mechanism
c = 0.95
# b: branching factor of the semantic version control
b = 3

# Model the FAIR components based on the parameters
# Findability (f) and Accessibility (a) are primarily limited by the consistency level.
f = c
a = c

# Interoperability (i) and Reusability (r) are limited by both consistency
# and the complexity introduced by the semantic branching factor.
i = c / b
r = c / b

# Calculate the overall FAIR compliance score (R) as the arithmetic mean of the components.
R = (f + a + i + r) / 4

# Output the explanation and the final equation with calculated values.
print("The FAIR compliance score (R) is modeled as the average of its four components.")
print(f"Given c = {c} and b = {b}:")
print(f"f = c = {f}")
print(f"a = c = {a}")
print(f"i = c / b = {c} / {b} = {i}")
print(f"r = c / b = {c} / {b} = {r}")
print("\nThe final equation for R is (f + a + i + r) / 4:")
print(f"R = ({f} + {a} + {i} + {r}) / 4 = {R}")
