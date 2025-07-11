import sys

# Step 1: Define the function for the Horizontal Transformation
def horizontal_transform(triplet):
    """Applies the horizontal transformation rules to a triplet."""
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

# Step 2: Define the known starting triplet for the third row
triplet_3_1 = [7, 2, 9]

# Step 3: Calculate the first missing triplet (the middle one in the last row)
print(f"Calculating the first missing triplet from {triplet_3_1}:")
missing_triplet_1 = horizontal_transform(triplet_3_1)
print(f"The first missing triplet is: {missing_triplet_1}\n")

# Step 4: Calculate the second missing triplet (the rightmost one)
print(f"Calculating the second missing triplet from {missing_triplet_1}:")
missing_triplet_2 = horizontal_transform(missing_triplet_1)
print(f"The second missing triplet is: {missing_triplet_2}\n")

# Step 5: Collect all missing elements and calculate their sum
missing_elements = missing_triplet_1 + missing_triplet_2
total_sum = sum(missing_elements)

# Step 6: Print the final summation equation as requested
equation_parts = [str(num) for num in missing_elements]
equation_str = " + ".join(equation_parts)

print("The final sum of the missing elements is calculated as follows:")
print(f"{equation_str} = {total_sum}")

# The final numeric answer
sys.stdout.write(f"\n<<<24>>>")