import math

# --- Given Matrix Data ---
L2 = [5, 7, 2]
M2 = [8, 4, 10]
R2 = [3, 1, 8]
L3 = [7, 2, 9]

# --- Helper function to check for prime numbers ---
def is_prime(n):
    return n in [2, 3, 5, 7, 11]

# --- Transformation Rules Implementation ---

def horizontal_transform(triplet):
    """Applies the horizontal transformation rule."""
    x, y, z = triplet
    if x + y > 10:
        next_x = (x * 3 - y) % 12
        next_y = (y * 2 + 4) % 12
        next_z = (z + x) % 12
    else:
        next_x = (x * 2 + y) % 12
        next_y = (y * 3 - 2) % 12
        next_z = (z * 2) % 12
    return [next_x, next_y, next_z]

def vertical_transform(temp_triplet, prev_triplet):
    """Applies the vertical transformation rule."""
    x, y, z = temp_triplet
    prev_x, prev_y, prev_z = prev_triplet
    if is_prime(prev_z):
        next_x = (x - 3 + prev_y) % 12
        next_y = (y + prev_x) % 12
        next_z = (z * 2 + prev_x) % 12
    else:
        next_x = (x + 2 - prev_y) % 12
        next_y = (y * 2 - prev_x) % 12
        next_z = (z + 3 + prev_z) % 12
    return [next_x, next_y, next_z]

# --- Main Calculation Logic ---

# 1. Calculate the middle triplet of row 3 (M3)
# Horizontal step using L3
temp_M3 = horizontal_transform(L3)
# Vertical step using the temp triplet and M2
M3 = vertical_transform(temp_M3, M2)
print(f"Calculated Middle Triplet of Row 3: {M3}")

# 2. Calculate the right triplet of row 3 (R3)
# Horizontal step using the newly calculated M3
temp_R3 = horizontal_transform(M3)
# Vertical step using the temp triplet and R2
R3 = vertical_transform(temp_R3, R2)
print(f"Calculated Right Triplet of Row 3: {R3}")

# 3. Sum the missing elements
missing_elements = M3 + R3
total_sum = sum(missing_elements)

# 4. Print the final result in the specified format
equation_parts = [str(n) for n in missing_elements]
equation_str = " + ".join(equation_parts)
print(f"\nThe sum of the missing elements is calculated as follows:")
print(f"{equation_str} = {total_sum}")

print(f"\n<<<The sum of the missing elements is {total_sum}>>>")
