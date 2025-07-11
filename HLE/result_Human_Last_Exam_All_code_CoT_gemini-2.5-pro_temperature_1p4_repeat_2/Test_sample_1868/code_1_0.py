# Given parameters for the federated knowledge graph system
# c: consistency level of the decentralized identifier resolution
c = 0.95
# b: branching factor of semantic version control
b = 3

# --- Step 1: Model the FAIR components based on the parameters ---

# Findability (f) is limited by the identifier consistency level.
# Maximum theoretical findability is assumed to be equal to c.
f = c

# Accessibility (a) is assumed to be perfect under the "best implementation" clause.
# Maximum theoretical accessibility is 1.
a = 1.0

# Interoperability (i) is negatively affected by the branching factor.
# We model this as an inverse relationship, i = 1/b.
i = 1 / b

# Reusability (r) is also assumed to be perfect.
# Maximum theoretical reusability is 1.
r = 1.0

# --- Step 2: Calculate the overall FAIR score (R) ---

# The overall FAIR score R is the arithmetic mean of the four components.
R = (f + a + i + r) / 4

# --- Step 3: Print the final equation and result ---

# We print the calculation step-by-step to show how the final score is derived.
print("Calculating the theoretical maximum FAIR compliance score (R):")
print(f"R = (f + a + i + r) / 4")
# The f-string formats the floating point numbers to two decimal places for clarity in the equation.
print(f"R = ({f:.2f} + {a:.2f} + {i:.2f} + {r:.2f}) / 4")
print(f"R = {f + a + i + r:.3f} / 4")
print(f"The final score is R = {R:.4f}")

# The final answer in the required format
print(f"\n<<<{R:.4f}>>>")