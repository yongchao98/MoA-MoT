def solve_ratio():
    """
    This function solves for the ratio BM/MI in terms of the side lengths a, b, c
    of triangle ABC and demonstrates the result with a numerical example.

    The derivation is as follows:
    1.  Let the side lengths opposite to vertices A, B, C be a, b, c.
    2.  The point M, where the angle bisector of angle B intersects the circumcircle,
        is equidistant from A, C, and the incenter I. This gives the key equality: MA = MC = MI.
    3.  Applying Ptolemy's Theorem to the cyclic quadrilateral ABCM, we get:
        AC * BM = AB * CM + BC * AM
    4.  Substituting the known lengths:
        b * BM = c * CM + a * AM
    5.  Using the equality from step 2 (CM = AM = MI):
        b * BM = c * MI + a * MI
        b * BM = (a + c) * MI
    6.  Rearranging the terms gives the final ratio:
        BM / MI = (a + c) / b
    """
    
    # We can demonstrate this formula with a sample triangle.
    # Let's use a triangle with side lengths a = 7, b = 8, c = 9.
    # These values satisfy the triangle inequality (e.g., 7+8 > 9).
    
    a = 7
    b = 8
    c = 9
    
    # Calculate the ratio using the derived formula
    ratio = (a + c) / b
    
    # Output the results in a clear format.
    print("The derived symbolic formula for the ratio BM/MI is:")
    print("BM / MI = (a + c) / b")
    
    print("\n--- Example Calculation ---")
    print(f"For a triangle with side lengths a = {a}, b = {b}, c = {c}:")
    
    # Printing the final equation with the numbers plugged in.
    print("\nThe final equation and calculated ratio are:")
    # The string representation shows the formula with the specific numbers
    equation_str = f"BM / MI = ({a} + {c}) / {b}"
    print(f"{equation_str} = {ratio}")

# Execute the function to see the output
solve_ratio()
