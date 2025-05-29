# Initial quantities
crop_A = 4
crop_B = 2
product_X = 0

# Function to apply the methods in cycle
def apply_methods(crop_A, crop_B, product_X):
    while True:
        # Try Method 1
        if crop_A >= 1 and crop_B >= 1:
            crop_A -= 1
            crop_B -= 1
            product_X += 1
        else:
            # If Method 1 cannot be applied, try Method 2
            if crop_A >= 3:
                crop_A -= 3
                product_X += 2
            else:
                # If neither method can be applied, break the loop
                break
    return crop_A, crop_B, product_X

# Apply the methods
remaining_A, remaining_B, obtained_X = apply_methods(crop_A, crop_B, product_X)

# Output the final result
print([remaining_A, remaining_B, obtained_X])