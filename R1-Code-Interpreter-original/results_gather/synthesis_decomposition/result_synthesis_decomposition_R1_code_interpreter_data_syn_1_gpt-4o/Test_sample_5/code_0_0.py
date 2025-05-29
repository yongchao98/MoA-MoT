# Initial quantities
crop_A = 5
crop_B = 3
product_X = 0

# Function to execute the methods
def exchange_crops(crop_A, crop_B, product_X):
    while True:
        # Method 1: 1 unit of A and 1 unit of B for 1 unit of X
        if crop_A >= 1 and crop_B >= 1:
            crop_A -= 1
            crop_B -= 1
            product_X += 1
        else:
            # If Method 1 cannot be executed, try Method 2
            # Method 2: 3 units of A for 2 units of X
            if crop_A >= 3:
                crop_A -= 3
                product_X += 2
            else:
                # If neither method can be executed, break the loop
                break
    return crop_A, crop_B, product_X

# Execute the exchange process
remaining_A, remaining_B, obtained_X = exchange_crops(crop_A, crop_B, product_X)

# Output the final result
print([remaining_A, remaining_B, obtained_X])