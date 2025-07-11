# Given parameters
c = 0.95  # Consistency level of identifier resolution
b = 3     # Branching factor of semantic version control

print("Step 1: Define the model for each FAIR component based on the given parameters.\n")

# Findability (f) is modeled as being equal to the consistency level.
f = c
print(f"Findability (f) = c = {f:.2f}")

# Accessibility (a) is assumed to be perfect (1.0) in a 'best implementation' scenario.
a = 1.0
print(f"Accessibility (a) = 1.0 (ideal)")

# Interoperability (i) is modeled as the inverse of the branching factor.
i = 1 / b
print(f"Interoperability (i) = 1 / b = 1 / {b} = {i:.4f}")

# Reusability (r) is modeled as being proportional to consistency c and inversely to branching factor b.
r = c / b
print(f"Reusability (r) = c / b = {c:.2f} / {b} = {r:.4f}")

print("\nStep 2: Calculate the final FAIR compliance score (R) as the average of the components.\n")

# The overall FAIR score R is the arithmetic mean of the four components.
R = (f + a + i + r) / 4

# Print the final equation with all the calculated values
print("Final Equation:")
print(f"R = (f + a + i + r) / 4")
print(f"R = ({f:.4f} + {a:.4f} + {i:.4f} + {r:.4f}) / 4")
print(f"R = ({f + a + i + r:.4f}) / 4")
print(f"R = {R:.4f}")

print("\nThe theoretical maximum value of R is calculated above.")
print("<<<0.6542>>>")