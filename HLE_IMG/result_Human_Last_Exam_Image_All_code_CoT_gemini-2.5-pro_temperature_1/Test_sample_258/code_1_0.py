import math

def calculate_cp_length():
    """
    Calculates the length of segment CP based on the given length of AB.
    
    The geometric construction leads to the relationship:
    CP = AB * (1 + sqrt(5)) / 2
    This is the golden ratio multiplied by the length of AB.
    """
    
    # The given length of the segment AB.
    ab_length = 18468
    
    # The value of the square root of 5.
    sqrt_5 = math.sqrt(5)
    
    # Calculate the length of CP using the derived formula.
    cp_length = ab_length * (1 + sqrt_5) / 2
    
    # The problem requires outputting each number in the final equation.
    # The following print statement shows the full calculation with numerical values.
    print(f"CP = {ab_length} * (1 + {sqrt_5}) / 2 = {cp_length:.4f}")

if __name__ == "__main__":
    calculate_cp_length()