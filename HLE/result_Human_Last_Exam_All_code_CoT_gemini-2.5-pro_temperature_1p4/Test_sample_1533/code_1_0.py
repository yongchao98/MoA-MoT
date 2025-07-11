def solve_triangle_ratio():
    """
    This function provides the solution to the geometry problem by printing
    the derived formula for the ratio BM/MI in terms of the side lengths a, b, and c.
    """
    
    # In triangle ABC, the side lengths are denoted by a, b, and c, where:
    # a is the length of side BC (opposite to angle A)
    # b is the length of side CA (opposite to angle B)
    # c is the length of side AB (opposite to angle C)
    
    # We will represent these side lengths as symbolic strings for the output.
    side_a = "a"
    side_b = "b"
    side_c = "c"
    
    # Based on the detailed geometric derivation, the ratio BM/MI is found to be (a + c) / b.
    # We will now construct and print this expression as the final answer.
    
    numerator = f"{side_a} + {side_c}"
    denominator = f"{side_b}"
    
    print("The final expression for the ratio BM/MI in terms of the side lengths a, b, and c is:")
    
    # Printing the final equation with its components as requested.
    print(f"BM / MI = ({numerator}) / {denominator}")

# Execute the function to display the solution.
solve_triangle_ratio()