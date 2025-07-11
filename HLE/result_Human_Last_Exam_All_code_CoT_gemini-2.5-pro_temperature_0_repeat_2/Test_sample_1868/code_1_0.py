import sys

# Step 1: Define the given parameters.
# c: Consistency level of the decentralized identifier resolution mechanism.
c = 0.95
# b: Branching factor of semantic version control.
b = 3

# Step 2: Model the FAIR metrics based on the parameters and assumptions.
# f (Findability) is modeled as being limited by the consistency level 'c'.
f = c
# a (Accessibility) is assumed to be perfect (1.0) in a best-case scenario.
a = 1.0
# i (Interoperability) is modeled as inversely proportional to the branching factor 'b'.
i = 1 / b
# r (Reusability) is also modeled as inversely proportional to the branching factor 'b'.
r = 1 / b

# Step 3: Calculate the overall FAIR compliance score R as the arithmetic mean.
R = (f + a + i + r) / 4

# Step 4: Print the final equation with each number and the result.
# The sys.stdout.write is used to prevent adding an extra newline at the end.
print("The theoretical maximum value of R is calculated using the following model:")
print(f"R = (Findability + Accessibility + Interoperability + Reusability) / 4")
print(f"R = (c + 1.0 + 1/b + 1/b) / 4")
print(f"R = ({f} + {a} + 1/{b} + 1/{b}) / 4")
print(f"R = ({f} + {a} + {i:.4f} + {r:.4f}) / 4")
sys.stdout.write(f"R = {R:.4f}\n")