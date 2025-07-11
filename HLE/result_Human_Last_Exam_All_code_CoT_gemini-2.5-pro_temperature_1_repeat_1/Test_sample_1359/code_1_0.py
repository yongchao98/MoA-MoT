import math

def is_prime(n):
    """Checks if a number is prime. For this problem, numbers are 0-11."""
    if n <= 1:
        return False
    for i in range(2, int(math.sqrt(n)) + 1):
        if n % i == 0:
            return False
    return True

def apply_horizontal_transform(triplet):
    """Applies the horizontal transformation rules."""
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

def apply_vertical_transform(current_triplet, previous_triplet):
    """Applies the vertical transformation rules."""
    x, y, z = current_triplet
    prev_x, prev_y, prev_z = previous_triplet

    if is_prime(prev_z):
        next_x = (x - 3 + prev_y) % 12
        next_y = (y + prev_x) % 12
        next_z = (z * 2 + prev_x) % 12
    else:
        next_x = (x + 2 - prev_y) % 12
        next_y = (y * 2 - prev_x) % 12
        next_z = (z + 3 + prev_z) % 12
    return [next_x, next_y, next_z]

# Known triplets from the matrix
L3 = [7, 2, 9]  # Left, Row 3
M2 = [8, 4, 10] # Middle, Row 2
R2 = [3, 1, 8]  # Right, Row 2

# --- Step 1: Calculate the middle triplet of Row 3 ([? ? ?]) ---

# First, apply horizontal transformation from L3
M3_intermediate = apply_horizontal_transform(L3)

# Then, apply vertical transformation using M2
M3_final = apply_vertical_transform(M3_intermediate, M2)

# --- Step 2: Calculate the right triplet of Row 3 ([? ? ?]) ---

# First, apply horizontal transformation from the newly found M3_final
R3_intermediate = apply_horizontal_transform(M3_final)

# Then, apply vertical transformation using R2
R3_final = apply_vertical_transform(R3_intermediate, R2)

# --- Step 3: Sum the elements and print the result ---
missing_elements = M3_final + R3_final
total_sum = sum(missing_elements)

print(f"The first missing triplet is: {M3_final}")
print(f"The second missing triplet is: {R3_final}")
equation_str = " + ".join(map(str, missing_elements))
print(f"The sum of the missing elements is {equation_str} = {total_sum}")

<<<20>>>