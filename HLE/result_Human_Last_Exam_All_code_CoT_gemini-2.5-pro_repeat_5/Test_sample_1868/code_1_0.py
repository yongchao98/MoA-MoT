import math

# --- Step 1: Define the given parameters ---
# Consistency level of the decentralized identifier resolution mechanism
c = 0.95
# Branching factor of semantic version control
b = 3

# --- Step 2: Model the theoretical maximum for each FAIR metric ---

# Findability (f) is limited by identifier consistency.
f = c

# Accessibility (a) is assumed to be perfect (1.0) under best practices.
a = 1.0

# Interoperability (i) is inversely affected by the semantic branching factor.
i = 1 / b

# Reusability (r) is limited by provenance, which depends on identifier consistency.
r = c

# --- Step 3: Calculate the overall FAIR compliance score (R) ---
# R is the arithmetic mean of the four metrics.
R = (f + a + i + r) / 4

# --- Step 4: Output the equation and the final result ---
print("This script calculates the theoretical maximum FAIR compliance score (R).")
print("\n--- Model ---")
print(f"Findability (f) = Consistency (c) = {f}")
print(f"Accessibility (a) = Best Practice = {a}")
print(f"Interoperability (i) = 1 / Branching Factor (b) = {i:.4f}")
print(f"Reusability (r) = Consistency (c) = {r}")

print("\n--- Calculation ---")
print(f"R = (f + a + i + r) / 4")
print(f"R = ({f} + {a} + 1/{b} + {r}) / 4")
print(f"R = ({f} + {a} + {i:.4f} + {r}) / 4")
print(f"R = ({f + a + i + r:.4f}) / 4")

print("\n--- Result ---")
# Using 5 decimal places for a more precise final answer
print(f"The theoretical maximum value of R is: {R:.5f}")
