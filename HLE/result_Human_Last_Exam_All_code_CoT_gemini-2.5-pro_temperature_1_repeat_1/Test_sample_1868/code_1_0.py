import sys

# Step 1: Define the given parameters.
# c: Consistency level of the decentralized identifier resolution mechanism.
c = 0.95
# b: Branching factor of the semantic version control.
b = 3

# Step 2: Model the FAIR components based on the parameters.
# Findability (f) and Accessibility (a) are limited by the identifier consistency 'c'.
f = c
a = c

# Interoperability (i) and Reusability (r) are limited by the complexity
# introduced by the branching factor 'b'. We model this as an inverse relationship.
i = 1 / b
r = 1 / b

# Step 3: Calculate the overall FAIR compliance score R as the average of its components.
R = (f + a + i + r) / 4

# Step 4: Output the results as requested.
# The final equation is printed with each number substituted.
print("The theoretical maximum FAIR score (R) is calculated as the average of its four components.")
print("Model:")
print(f"  Findability (f) = c = {f}")
print(f"  Accessibility (a) = c = {a}")
print(f"  Interoperability (i) = 1/b = {i}")
print(f"  Reusability (r) = 1/b = {r}")
print("\nFinal Equation:")
# The following line prints the equation with the component values, as requested.
print(f"R = ({f} + {a} + {i} + {r}) / 4")

print(f"\nCalculated Value:")
print(f"R = {R}")

# Appending the final answer in the specified format to the output stream.
# This part will not be visible in standard execution but is for the final answer extraction.
sys.stdout.write(f"\n<<<{R}>>>")