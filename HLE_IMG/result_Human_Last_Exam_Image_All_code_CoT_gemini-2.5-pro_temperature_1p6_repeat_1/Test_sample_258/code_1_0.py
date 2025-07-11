import math

def calculate_cp_length(ab_length):
    """
    Calculates the length of CP based on a specific geometric construction.

    The derivation using coordinate geometry shows that the relationship between
    the length of CP and the length of AB is:
    CP = AB * (1 + sqrt(5)) / 2
    This is equivalent to AB * phi, where phi is the golden ratio.
    """
    
    # The given length of the initial segment AB.
    ab = float(ab_length)
    
    # Calculate the value of the square root of 5.
    sqrt_5 = math.sqrt(5)
    
    # The term (1 + sqrt(5)) / 2 is the golden ratio, phi.
    golden_ratio = (1 + sqrt_5) / 2
    
    # Calculate the final length of CP.
    cp_length = ab * golden_ratio
    
    print(f"Given the length of AB = {ab}")
    print("The final formula derived from the geometry is: CP = AB * (1 + sqrt(5)) / 2")
    print(f"Step 1: Calculate the value of (1 + sqrt(5)) / 2")
    print(f"   (1 + {sqrt_5}) / 2 = {golden_ratio}")
    print(f"Step 2: Calculate the length of CP")
    print(f"   CP = {ab} * {golden_ratio}")
    print(f"Final Answer: The length of CP is {cp_length:.4f}")

    return cp_length

# The length of AB provided in the problem.
AB_given_length = 18468

# Calculate and print the result.
result = calculate_cp_length(AB_given_length)

# Output the final answer in the required format.
print(f"<<<{result:.4f}>>>")
