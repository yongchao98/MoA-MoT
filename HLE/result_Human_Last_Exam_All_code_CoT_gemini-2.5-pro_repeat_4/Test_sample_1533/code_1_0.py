def solve_triangle_ratio():
    """
    This function calculates the ratio BM/MI for a triangle ABC with given side lengths.
    It demonstrates the calculation for a sample triangle and prints the general formula.
    
    The side lengths a, b, c are opposite to vertices A, B, C respectively.
    """
    
    # Example side lengths for a triangle (e.g., a 5-6-7 triangle)
    # To form a valid triangle, the sum of any two sides must be greater than the third.
    a = 5
    b = 6
    c = 7
    
    print(f"Consider a triangle ABC with side lengths a = {a}, b = {b}, c = {c}.")
    
    # The derived formula for the ratio BM/MI is (a + c) / b
    
    # Calculate the numerator and denominator of the ratio
    ratio_numerator = a + c
    ratio_denominator = b
    
    # Calculate the final value of the ratio
    if ratio_denominator == 0:
        print("Error: Side length b cannot be zero.")
        return
        
    ratio_value = ratio_numerator / ratio_denominator
    
    print("\nThe ratio BM/MI can be expressed in terms of the side lengths a, b, and c.")
    print("The derived formula is:")
    print("  BM     a + c")
    print("---- = -------")
    print("  MI       b")
    
    print("\nFor the given triangle:")
    print(f"  BM     {a} + {c}     {ratio_numerator}")
    print(f"---- = ------- = --- = {ratio_value:.4f}")
    print("  MI       {b}        {ratio_denominator}")

# Execute the function to see the output
solve_triangle_ratio()
