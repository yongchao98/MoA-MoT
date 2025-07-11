import math

def calculate_cp_length(ab_length):
    """
    Calculates the length of segment CP based on the length of segment AB.

    Args:
        ab_length (float): The length of the segment AB.

    Returns:
        float: The length of the segment CP.
    """
    # The geometric analysis shows that the length of CP is related to AB
    # by the golden ratio, phi.
    # The formula is: CP = AB * ( (sqrt(5) + 1) / 2 )
    
    # Calculate the components of the formula
    sqrt_5 = math.sqrt(5)
    phi = (sqrt_5 + 1) / 2
    cp_length = ab_length * phi
    
    # Print the steps of the calculation as requested
    print("The relationship between CP and AB is given by the golden ratio (phi):")
    print("Formula: CP = AB * ( (sqrt(5) + 1) / 2 )")
    print("\nGiven AB = {}".format(ab_length))
    print("We substitute the values into the formula:")
    
    # Showing each number in the final equation
    print("CP = {} * ( ({} + 1) / 2 )".format(ab_length, sqrt_5))
    print("CP = {} * ( {} / 2 )".format(ab_length, sqrt_5 + 1))
    print("CP = {} * {}".format(ab_length, phi))
    print("CP = {}".format(cp_length))
    
    # Print the final answer formatted to 4 decimal places
    print("\nThe final length of CP rounded to 4 decimal places is: {:.4f}".format(cp_length))
    
    return cp_length

# The length of AB is given as 18468
AB = 18468
final_cp_length = calculate_cp_length(AB)

# The problem asks for the answer in a specific format at the end.
# The following line is for the final answer extraction and is commented out
# to prevent duplicate output in a standard execution environment.
# print("<<< {:.4f} >>>".format(final_cp_length))