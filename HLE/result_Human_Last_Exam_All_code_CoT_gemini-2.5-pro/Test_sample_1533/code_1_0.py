def solve_ratio():
    """
    Calculates the ratio BM/MI for a triangle ABC with side lengths a, b, c.
    The formula derived is (a + c) / b.
    This function demonstrates the result with an example triangle.
    """
    # Example side lengths for a valid triangle (e.g., a Heron triangle)
    a = 7
    b = 8
    c = 9

    # The problem asks for the ratio BM/MI in terms of a, b, and c.
    # Based on the geometric derivation, the final expression is (a + c) / b.
    
    print("In a triangle ABC with side lengths a, b, and c:")
    print(f"Let side a = {a}")
    print(f"Let side b = {b}")
    print(f"Let side c = {c}")
    print("-" * 20)
    
    print("The ratio BM/MI is given by the formula: (a + c) / b")
    
    # Calculate the values for the formula
    numerator = a + c
    denominator = b
    ratio = numerator / denominator
    
    # Print the step-by-step calculation with the given numbers
    print("\nCalculation steps for the example:")
    print(f"BM / MI = ({a} + {c}) / {b}")
    print(f"BM / MI = {numerator} / {denominator}")
    print(f"BM / MI = {ratio}")

# Execute the function to see the output
solve_ratio()
