import sys

# Given parameters
# c: consistency level of the decentralized identifier resolution mechanism
c = 0.95
# b: branching factor of semantic version control
b = 3

# --- Step-by-step derivation of the FAIR score R ---

# The FAIR compliance score R is defined as the average of its four components:
# Findability (f), Accessibility (a), Interoperability (i), and Reusability (r).
# R = (f + a + i + r) / 4

# We calculate the theoretical maximum for each component based on the given parameters
# and the assumption of "best implementation of FAIR practices".

# 1. Maximum Findability (f_max)
# Findability is fundamentally limited by the consistency 'c' of the identifier
# resolution system. A perfect score of 1.0 is not possible if identifiers
# do not always resolve correctly.
f_max = c

# 2. Maximum Accessibility (a_max)
# Assuming best practices, accessibility protocols are open and standardized.
# The given parameters do not impose a theoretical limit on accessibility.
a_max = 1.0

# 3. Maximum Interoperability (i_max)
# Interoperability is hindered by semantic divergence. A branching factor 'b' > 1
# means concepts can have multiple versions, making integration difficult.
# The maximum interoperability is inversely proportional to this branching factor.
i_max = 1 / b

# 4. Maximum Reusability (r_max)
# Reusability requires clear provenance and metadata, which depends on both
# reliable identifiers (c) and a clear version history (b). It is therefore
# limited by both factors. A simple model is the ratio of c to b.
r_max = c / b

# Calculate the final theoretical maximum FAIR score (R_max)
R_max = (f_max + a_max + i_max + r_max) / 4

# --- Output the results step-by-step ---

print("Calculating the theoretical maximum FAIR compliance score (R_max)")
print("="*60)
print(f"Given parameters:\n  - Identifier Consistency (c) = {c}\n  - Branching Factor (b) = {b}\n")

print("The overall score R is the average of its four components:")
print("R = (Findability + Accessibility + Interoperability + Reusability) / 4\n")

print("Step 1: Calculate maximum Findability (f_max)")
print(f"  f_max = c = {f_max}\n")

print("Step 2: Determine maximum Accessibility (a_max)")
print(f"  a_max = 1.0 (based on best practices assumption)\n")

print("Step 3: Calculate maximum Interoperability (i_max)")
print(f"  i_max = 1 / b = 1 / {b} = {i_max:.4f}\n")

print("Step 4: Calculate maximum Reusability (r_max)")
print(f"  r_max = c / b = {c} / {b} = {r_max:.4f}\n")

print("Step 5: Substitute values into the final equation for R_max")
print(f"  R_max = (f_max + a_max + i_max + r_max) / 4")
print(f"  R_max = ({f_max} + {a_max} + {i_max:.4f} + {r_max:.4f}) / 4")
sum_of_components = f_max + a_max + i_max + r_max
print(f"  R_max = {sum_of_components:.4f} / 4")
print("-" * 60)
print(f"Final Theoretical Maximum Score (R_max) = {R_max:.4f}")
print("="*60)

# The following line is for programmatic parsing of the final answer.
# It should be the last line of the output.
sys.stdout.write(f"\n<<<{R_max:.2f}>>>")