import sys

# Define the input parameters
# c: Consistency level of the decentralized identifier resolution
c = 0.95
# b: Branching factor of semantic version control
b = 3

# --- Model Formulation ---
# The overall FAIR score (R) is modeled as the average of the four components.
# We assume a "best implementation," meaning each component score is 1.0
# unless limited by the given parameters c and b.

# 1. Findability (f)
# Findability is directly limited by the consistency of identifier resolution.
# If an identifier doesn't resolve correctly, the data cannot be found.
f = c

# 2. Accessibility (a)
# No specific parameters limit accessibility. With best practices (e.g., standard
# protocols, open access), we assume it is ideal.
a = 1.0

# 3. Interoperability (i)
# Interoperability depends on both resolving links to other data (c) and
# semantic compatibility. With 'b' branches, the chance of two systems using
# a compatible semantic version is modeled as 1/b.
i = c * (1 / b)

# 4. Reusability (r)
# Reusability is diminished by semantic ambiguity and complex provenance.
# A higher branching factor 'b' makes provenance tracking harder and introduces
# ambiguity, thus reducing reusability. We model this impact as 1/b.
r = 1 / b

# Calculate the final FAIR score R
R = (f + a + i + r) / 4

# --- Output the results ---
print("Calculating the theoretical maximum FAIR compliance score (R)...")
print("-" * 60)
print(f"Given Parameters:")
print(f"  - Identifier Consistency (c): {c}")
print(f"  - Branching Factor (b): {b}\n")

print("Derived FAIR Component Scores:")
print(f"  1. Findability (f) = c                                 = {f:.3f}")
print(f"  2. Accessibility (a) = Ideal                           = {a:.3f}")
print(f"  3. Interoperability (i) = c * (1/b)                    = {i:.3f}")
print(f"  4. Reusability (r) = 1/b                             = {r:.3f}\n")

print("Final Score Calculation (R):")
print("R = (f + a + i + r) / 4")
# Using format to explicitly show each number in the final equation
print(f"R = ({f:.3f} + {a:.3f} + {i:.3f} + {r:.3f}) / 4")
print(f"R = {f+a+i+r:.3f} / 4")
print(f"R = {R:.3f}")
print("-" * 60)

# The following line is used for programmatic retrieval of the answer.
# It prints the final numerical answer to stdout.
sys.stdout.write(f"\n<<<__{R:.2f}__>>>\n")