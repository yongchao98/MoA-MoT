import math

# Define the given parameters
# c: Consistency level of the decentralized identifier resolution mechanism
c = 0.95
# b: Branching factor of the semantic version control
b = 3

# --- Model the FAIR components based on the plan ---

# f (Findability) is modeled as being limited by the consistency 'c'.
f = c

# a (Accessibility) is assumed to be ideal (1.0) under the "best implementation".
a = 1.0

# i (Interoperability) is modeled as inversely proportional to the branching factor 'b'.
i = 1 / b

# r (Reusability) is also modeled as inversely proportional to the branching factor 'b'.
r = 1 / b

# --- Calculate the theoretical maximum FAIR score R ---

# R is the sum of the individual component scores.
R = f + a + i + r

# --- Output the results ---
print("This script calculates the theoretical maximum FAIR compliance score (R).")
print("The model assumes R is the sum of scores for Findability(f), Accessibility(a), Interoperability(i), and Reusability(r).")
print("-" * 30)
print(f"Given parameters:")
print(f"Identifier Consistency (c) = {c}")
print(f"Branching Factor (b) = {b}")
print("-" * 30)
print("Calculating component scores:")
print(f"f = c = {f:.2f}")
print(f"a = Best Practice Assumption = {a:.2f}")
print(f"i = 1/b = 1/{b} = {i:.4f}")
print(f"r = 1/b = 1/{b} = {r:.4f}")
print("-" * 30)
print("Final Equation:")
print(f"R = f + a + i + r")
print(f"R = {f:.2f} + {a:.2f} + {i:.4f} + {r:.4f}")
print(f"R = {R:.4f}")