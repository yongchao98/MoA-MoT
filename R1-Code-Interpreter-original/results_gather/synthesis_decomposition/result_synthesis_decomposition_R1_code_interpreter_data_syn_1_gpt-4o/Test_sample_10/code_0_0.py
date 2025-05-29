# Initial quantities
crop_A = 4
crop_B = 3
product_X = 0

# Function to apply the methods
def apply_methods(crop_A, crop_B, product_X):
    while True:
        # Try Method 1
        if crop_A >= 1 and crop_B >= 1:
            crop_A -= 1
            crop_B -= 1
            product_X += 1
        else:
            break  # Cannot apply Method 1

        # Try Method 2
        if crop_A >= 3:
            crop_A -= 3
            product_X += 2
        else:
            break  # Cannot apply Method 2

    return crop_A, crop_B, product_X

# Apply the methods
remaining_A, remaining_B, obtained_X = apply_methods(crop_A, crop_B, product_X)

# Output the result
print([remaining_A, remaining_B, obtained_X])