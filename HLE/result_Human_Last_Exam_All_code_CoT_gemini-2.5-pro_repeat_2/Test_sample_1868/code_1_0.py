import math

# Given parameters
c = 0.95  # Consistency level of the identifier resolution mechanism
b = 3     # Branching factor of semantic version control

# --- Model Formulation ---
# Based on the plan, we define the FAIR components in terms of c and b.
# We assume the best implementation of FAIR practices, where the score
# for each component is at its theoretical maximum given the system's constraints.

# Findability (f) is limited by identifier resolution consistency.
f = c

# Accessibility (a) is also limited by identifier resolution consistency.
a = c

# Interoperability (i) is negatively impacted by the semantic divergence
# caused by the branching factor.
i = 1 / b

# Reusability (r) is dependent on finding the correct version (c) but
# is diluted by the number of branches (b).
r = c / b

# The final FAIR compliance score (R) is the scaled average of the components.
# We scale it by 10 to get a score out of 10.
R = 10 * (f + a + i + r) / 4

# --- Output the results ---
print("--- Calculating the FAIR Compliance Score (R) ---")
print(f"Given parameters:")
print(f"  Identifier Resolution Consistency (c): {c}")
print(f"  Semantic Versioning Branching Factor (b): {b}\n")

print("Calculated FAIR component scores:")
print(f"  Findability (f) = c = {f:.4f}")
print(f"  Accessibility (a) = c = {a:.4f}")
print(f"  Interoperability (i) = 1 / b = {i:.4f}")
print(f"  Reusability (r) = c / b = {r:.4f}\n")

print("The final equation for the FAIR score R is:")
print(f"R = 10 * (f + a + i + r) / 4")
print(f"R = 10 * ({f:.2f} + {a:.2f} + {i:.4f} + {r:.4f}) / 4")
print(f"R = 10 * ({f + a + i + r:.4f}) / 4")
print(f"R = {10 * (f + a + i + r):.4f} / 4\n")

print(f"The theoretical maximum value of R is: {R:.4f}")
print("<<<6.375>>>")