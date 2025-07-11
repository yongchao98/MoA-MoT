def calculate_bm_mi_ratio(a, b, c):
    """
    Calculates the ratio BM/MI for a triangle ABC with side lengths a, b, c.

    Args:
        a (float): The length of side BC.
        b (float): The length of side AC.
        c (float): The length of side AB.
    
    Returns:
        float: The calculated ratio BM/MI.
    """
    # The derivation shows the ratio BM/MI is equal to (a + c) / b.
    if b == 0:
        print("Error: Side b cannot be zero.")
        return None
        
    numerator = a + c
    denominator = b
    ratio = numerator / denominator
    
    print(f"Given a triangle with side lengths a = {a}, b = {b}, c = {c}:")
    print(f"The formula for the ratio BM / MI is (a + c) / b.")
    print(f"The numerator is a + c = {a} + {c} = {numerator}")
    print(f"The denominator is b = {denominator}")
    print(f"The final ratio is {numerator} / {denominator} = {ratio}")
    
    return ratio

if __name__ == '__main__':
    # Example using a valid triangle (e.g., sides 13, 15, 14)
    # Side a is opposite to vertex A, b to B, c to C.
    side_a = 13
    side_b = 15
    side_c = 14
    
    # Check if the sides can form a triangle (Triangle Inequality Theorem)
    if (side_a + side_b > side_c) and \
       (side_a + side_c > side_b) and \
       (side_b + side_c > side_a):
        calculate_bm_mi_ratio(side_a, side_b, side_c)
    else:
        print(f"The sides a={side_a}, b={side_b}, c={side_c} cannot form a triangle.")
