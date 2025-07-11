import math

# Step 1 & 2: Initial values and identifying the core issue.
# Total records are 720. Lost are 1/6, so readable are 720 * 5/6 = 600.
# The problem states "One-fourth use both 'Augustus' and 'Caesar'", implying b = 600/4 = 150.
# However, the constraint `sqrt(b) + 3*a = d` requires `sqrt(b)` to likely be an integer for 'a' and 'd' to be integer record counts.
# sqrt(150) is not an integer, indicating a conflict.

# Step 3: Hypothesis - the fractions are approximations, and b must be a perfect square.
# The closest perfect square to 150 is 144.
b = 144

# Step 4: Solve for a, c, and d.
# The constraint "single variants (c) form a perfect square whose root equals the number of lost records (a) divided by 4"
# gives us sqrt(c) = a / 4. This means a must be a multiple of 4.
# The fractional description for 'c' is 2/5 of records not using both names.
# Records not using both names = 600 - b = 600 - 144 = 456.
# Approximate c = (2/5) * 456 = 182.4.
# We need to find 'a' (a multiple of 4) where c = (a/4)^2 is close to 182.4.
# a^2 / 16 should be near 182.4 -> a^2 should be near 2918.4 -> a should be near sqrt(2918.4) = 54.02.
# The closest multiple of 4 to 54.02 is 52.
a = 52

# Now calculate c based on a=52
# sqrt(c) = 52 / 4 = 13
c = 13**2

# Now calculate d using the main constraint equation: sqrt(b) + 3a = d
d = math.sqrt(b) + 3 * a

# Convert d to an integer as it represents a number of records
d = int(d)

# Step 5: Determine the divisor 'e'
# 'e' is the number of distinct ways Augustus is named.
# 1. Dual naming ("Augustus" and "Caesar") - variable b
# 2. Single-variant ("Octavianus" or "Augustus") - variable c
# 3. Full imperial title ("Imperator Caesar Augustus") - variable d
# This gives 3 distinct patterns.
e = 3

# Step 6: Final calculation
product = a * b * c * d
division_result = product / e
final_answer = math.floor(division_result / 1000)

print(f"Found values are:")
print(f"a) Lost records (Octavius ones): {a}")
print(f"b) Documents with dual naming: {b}")
print(f"c) Single-variant documents: {c}")
print(f"d) Full imperial title documents: {d}")
print(f"The number of distinct naming patterns for Augustus is: {e}")

print(f"\nCalculation:")
print(f"({a} * {b} * {c} * {d}) / {e} / 1000")
print(f"= {product} / {e} / 1000")
print(f"= {division_result} / 1000")
print(f"= {division_result / 1000}")
print(f"The floor of the result is: {final_answer}")
