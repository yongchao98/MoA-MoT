# Initial quantities
crop_A = 5
crop_B = 2
product_X = 0

# Method application cycle
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

# Output the final state
print([crop_A, crop_B, product_X])