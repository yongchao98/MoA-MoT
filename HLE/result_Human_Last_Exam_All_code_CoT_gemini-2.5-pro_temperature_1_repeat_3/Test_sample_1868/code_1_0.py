import math

# Given parameters
# c: Consistency level of the decentralized identifier resolution
c = 0.95
# b: Branching factor of semantic version control
b = 3

# Step 1: Model the maximum score for each FAIR component
# Findability (f) and Accessibility (a) are limited by the identifier consistency.
f_max = c
a_max = c

# Interoperability (i) and Reusability (r) are limited by semantic ambiguity
# from the branching factor. We model this as an inverse relationship.
i_max = 1 / b
r_max = 1 / b

# Step 2: Calculate the theoretical maximum FAIR compliance score (R_max)
# R is defined as the average of its four components.
R_max = (f_max + a_max + i_max + r_max) / 4

# Step 3: Print the final equation and the result
print("The theoretical maximum FAIR score (R_max) is calculated as the average of the maximum possible scores for each component.")
print("R_max = (f_max + a_max + i_max + r_max) / 4\n")
print("With the given constraints:")
print(f"f_max = c = {f_max}")
print(f"a_max = c = {a_max}")
print(f"i_max = 1/b = 1/{b} = {i_max:.4f}")
print(f"r_max = 1/b = 1/{b} = {r_max:.4f}\n")

print("Plugging the values into the equation:")
# Using format specifiers to show the numbers in the equation clearly
print(f"R_max = ({f_max} + {a_max} + {i_max:.4f} + {r_max:.4f}) / 4")
print(f"R_max = {f_max + a_max + i_max + r_max:.4f} / 4")
print(f"Final R_max = {R_max:.4f}")
<<<0.6417>>>