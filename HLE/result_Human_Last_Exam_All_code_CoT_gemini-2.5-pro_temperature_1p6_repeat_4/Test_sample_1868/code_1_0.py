import sys

# Given parameters
c = 0.95  # Consistency level of identifier resolution
b = 3     # Branching factor of semantic version control

# --- Step 1: Model the FAIR Components ---
# Findability (f) is modeled as the identifier consistency level.
f = c

# Accessibility (a) is assumed to be perfect under best practices.
a = 1.0

# Interoperability (i) is modeled as inversely proportional to the branching complexity.
# A higher branching factor makes interoperability harder.
i = 1 / b

# Reusability (r) is dependent on identifier consistency (provenance)
# but hindered by branching complexity.
r = c / b

# --- Step 2: Model the Overall FAIR Score (R) ---
# R is calculated as the average of the four components.
R = (f + a + i + r) / 4

# --- Step 3: Print the Calculation ---
# The request requires printing each number in the final equation.
print("Calculating the FAIR Compliance Score (R)")
print("-----------------------------------------")
print(f"Given inputs:")
print(f"  Identifier consistency (c): {c}")
print(f"  Branching factor (b): {b}")
print("\nModeling FAIR components:")
print(f"  Findability (f) = c = {f}")
print(f"  Accessibility (a) = 1.0 (assumed best practice)")
print(f"  Interoperability (i) = 1 / b = 1 / {b} = {i:.4f}")
print(f"  Reusability (r) = c / b = {c} / {b} = {r:.4f}")
print("\nFinal Equation for R = (f + a + i + r) / 4:")
# The following line explicitly prints each number in the equation.
print(f"R = ({f} + {a} + {round(i, 4)} + {round(r, 4)}) / 4 = {R:.2f}")

# To satisfy the format requirement, printing the final number alone.
# This part is for programmatic use if needed, but the answer is clear above.
# print(f"{R:.2f}")

sys.stdout.flush()
print("<<<0.65>>>")