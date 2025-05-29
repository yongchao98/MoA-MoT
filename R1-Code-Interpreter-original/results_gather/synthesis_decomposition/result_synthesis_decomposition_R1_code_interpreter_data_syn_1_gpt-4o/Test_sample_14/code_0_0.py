# Initial quantities
crop_A = 5
crop_B = 3
product_X = 0

# Simulate the process
while True:
    # Method 1
    if crop_A >= 1 and crop_B >= 1:
        crop_A -= 1
        crop_B -= 1
        product_X += 1
    else:
        break
    
    # Method 2
    if crop_A >= 3:
        crop_A -= 3
        product_X += 2
    else:
        break

# Output the final state
print([crop_A, crop_B, product_X])