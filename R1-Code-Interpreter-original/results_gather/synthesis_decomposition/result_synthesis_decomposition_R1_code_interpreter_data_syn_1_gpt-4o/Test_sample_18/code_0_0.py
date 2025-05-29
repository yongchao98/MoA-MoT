# Initial quantities
crop_A = 3
crop_B = 3
product_X = 0

# Method 1: 1A + 1B -> 1X
# Method 2: 3A -> 2X

while True:
    # Try Method 1
    if crop_A >= 1 and crop_B >= 1:
        crop_A -= 1
        crop_B -= 1
        product_X += 1
    else:
        break  # If Method 1 cannot be applied, break the loop

    # Try Method 2
    if crop_A >= 3:
        crop_A -= 3
        product_X += 2
    else:
        continue  # If Method 2 cannot be applied, continue to next cycle

# Output the final state
print([crop_A, crop_B, product_X])