def calculate_triangle_ratio(a, b, c):
    """
    Calculates the ratio BM/MI for a triangle ABC with side lengths a, b, c.

    The derivation shows that the ratio BM/MI is equal to (a + c) / b.
    This function takes the side lengths and computes this value.

    Args:
        a (float): The length of side BC, opposite to vertex A.
        b (float): The length of side AC, opposite to vertex B.
        c (float): The length of side AB, opposite to vertex C.
    """
    # Check if the side lengths can form a triangle
    if not (a + b > c and a + c > b and b + c > a):
        print(f"Error: The side lengths a={a}, b={b}, c={c} do not form a valid triangle.")
        return

    # The ratio BM/MI is derived to be (a + c) / b
    ratio = (a + c) / b

    print(f"For a triangle with side lengths a = {a}, b = {b}, c = {c}:")
    
    # The problem asks to output each number in the final equation.
    # The final equation is: Ratio = (a + c) / b
    print("The final equation for the ratio BM/MI is (a + c) / b.")
    print(f"Substituting the given values: ({a} + {c}) / {b}")
    print(f"The calculated ratio is: {ratio}")

# Example usage with a sample triangle (e.g., a=7, b=8, c=9)
# In this case, side 'a' is opposite to angle A, 'b' to B, and 'c' to C.
side_a = 7
side_b = 8
side_c = 9
calculate_triangle_ratio(side_a, side_b, side_c)