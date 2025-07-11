def is_prime(n):
    """Checks if a number is prime."""
    if n < 2:
        return False
    for i in range(2, int(n**0.5) + 1):
        if n % i == 0:
            return False
    return True

def calculate_next_triplet(prev_triplet):
    """
    Calculates the next triplet in a column using vertical transformation rules.
    """
    x_prev, y_prev, z_prev = prev_triplet
    
    if is_prime(z_prev):
        # Rule: If previous z is prime
        next_x = (x_prev - 3 + y_prev) % 12
        next_y = (y_prev + x_prev) % 12
        next_z = (z_prev * 2 + x_prev) % 12
    else:
        # Rule: If previous z is not prime
        next_x = (x_prev + 2 - y_prev) % 12
        next_y = (y_prev * 2 - x_prev) % 12
        next_z = (z_prev + 3 + z_prev) % 12
        
    return [next_x, next_y, next_z]

# Given triplets from Row 2
r2_middle = [8, 4, 10]
r2_right = [3, 1, 8]

# --- Step 1: Calculate Row 3, Middle Triplet ---
r3_middle = calculate_next_triplet(r2_middle)
x1, y1, z1 = r3_middle
print(f"Calculating the middle triplet of the third row from [8 4 10]:")
print(f"Previous z (10) is not prime.")
print(f"New x = (8 + 2 - 4) mod 12 = {x1}")
print(f"New y = (4 * 2 - 8) mod 12 = {y1}")
print(f"New z = (10 + 3 + 10) mod 12 = {z1}")
print(f"Resulting Triplet: [{x1} {y1} {z1}]")
print("-" * 20)

# --- Step 2: Calculate Row 3, Right Triplet ---
r3_right = calculate_next_triplet(r2_right)
x2, y2, z2 = r3_right
print(f"Calculating the right triplet of the third row from [3 1 8]:")
print(f"Previous z (8) is not prime.")
print(f"New x = (3 + 2 - 1) mod 12 = {x2}")
print(f"New y = (1 * 2 - 3) mod 12 = {y2}")
print(f"New z = (8 + 3 + 8) mod 12 = {z2}")
print(f"Resulting Triplet: [{x2} {y2} {z2}]")
print("-" * 20)

# --- Step 3: Calculate the sum of missing elements ---
missing_elements = r3_middle + r3_right
total_sum = sum(missing_elements)

print("The missing elements are:", ", ".join(map(str, missing_elements)))
equation = " + ".join(map(str, missing_elements))
print(f"The sum is calculated as: {equation} = {total_sum}")
