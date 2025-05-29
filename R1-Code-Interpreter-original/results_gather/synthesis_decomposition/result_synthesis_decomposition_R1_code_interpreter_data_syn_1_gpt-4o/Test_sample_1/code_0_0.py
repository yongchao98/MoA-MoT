# Initial quantities
crop_A = 4
crop_B = 2
product_X = 0

# Function to execute the methods
def execute_methods(crop_A, crop_B, product_X):
    # Method 1
    if crop_A >= 1 and crop_B >= 1:
        crop_A -= 1
        crop_B -= 1
        product_X += 1
    
    # Method 2
    if crop_A >= 3:
        crop_A -= 3
        product_X += 2
    
    return crop_A, crop_B, product_X

# Simulate the process
while True:
    new_A, new_B, new_X = execute_methods(crop_A, crop_B, product_X)
    # Check if any changes occurred
    if new_A == crop_A and new_B == crop_B and new_X == product_X:
        break
    crop_A, crop_B, product_X = new_A, new_B, new_X

# Output the final result
print([crop_A, crop_B, product_X])