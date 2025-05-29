# Initial quantities
crop_A = 3
crop_B = 3
product_X = 0

# Method execution loop
while True:
    # Try Method 1
    if crop_A >= 1 and crop_B >= 1:
        crop_A -= 1
        crop_B -= 1
        product_X += 1
    # Try Method 2
    elif crop_A >= 3:
        crop_A -= 3
        product_X += 2
    else:
        # If neither method can be executed, break the loop
        break

# Output the final state
print([crop_A, crop_B, product_X])