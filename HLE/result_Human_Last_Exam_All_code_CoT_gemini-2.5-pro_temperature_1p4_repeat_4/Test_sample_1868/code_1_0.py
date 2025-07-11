# Define the given system parameters
c = 0.95  # Identifier resolution consistency level
b = 3     # Semantic version control branching factor

# --- Step 1: Calculate the theoretical maximum score for each FAIR principle ---
# Findability and Accessibility scores are limited by the consistency 'c'.
# We assume a scale of 0-10.
f = 10 * c
a = 10 * c

# Interoperability and Reusability scores are limited by the branching factor 'b'.
# A higher 'b' complicates interoperability and reuse, so we model it as an inverse relationship.
i = 10 / b
r = 10 / b

print(f"Based on the model, the maximum sub-scores are:")
print(f"Findability (f) = 10 * c = {f:.2f}")
print(f"Accessibility (a) = 10 * c = {a:.2f}")
print(f"Interoperability (i) = 10 / b = {i:.2f}")
print(f"Reusability (r) = 10 / b = {r:.2f}")
print("-" * 30)

# --- Step 2: Calculate the final FAIR compliance score R ---
# R is the arithmetic mean of the four sub-scores.
R = (f + a + i + r) / 4

# --- Step 3: Output the final calculation and result ---
# Displaying the final equation with the calculated numbers
print("The final FAIR compliance score R is calculated as:")
print(f"R = (f + a + i + r) / 4")
print(f"R = ({f:.2f} + {a:.2f} + {i:.2f} + {r:.2f}) / 4")
total_score = f + a + i + r
print(f"R = {total_score:.2f} / 4")
print(f"The theoretical maximum value of R is: {R}")