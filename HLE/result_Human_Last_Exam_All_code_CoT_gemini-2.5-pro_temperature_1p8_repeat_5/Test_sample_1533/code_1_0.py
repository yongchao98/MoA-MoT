def calculate_bm_mi_ratio(a, b, c):
    """
    Calculates the ratio BM/MI for a triangle ABC with side lengths a, b, c.

    In triangle ABC, the angle bisectors of angle BAC and angle ABC intersect at point I.
    Let BI intersect the circumcircle of triangle ABC at a point M.
    This function calculates the ratio BM/MI.

    Args:
        a: Length of side BC (opposite to vertex A).
        b: Length of side AC (opposite to vertex B).
        c: Length of side AB (opposite to vertex C).
    
    Returns:
        The numerical value of the ratio BM/MI.
    """
    # The sides must form a valid triangle
    if not (a + b > c and a + c > b and b + c > a):
        print("Error: The given side lengths do not form a valid triangle.")
        return None

    # The ratio is given by the formula (a+c)/b
    ratio = (a + c) / b
    
    # Print the equation with the given numbers
    print(f"For a triangle with side lengths a={a}, b={b}, c={c}:")
    print("The ratio BM/MI is calculated using the formula (a + c) / b.")
    print(f"({a} + {c}) / {b} = {ratio}")
    return ratio

if __name__ == '__main__':
    # Example usage with a common triangle (e.g. 3, 5, 4 where b is hypotenuse).
    # In our problem context, a and c are the sides adjacent to angle B.
    # Let's take a=3, c=4, b=5. This makes angle B = 90 degrees.
    side_a = 3
    side_b = 5
    side_c = 4
    calculate_bm_mi_ratio(side_a, side_b, side_c)

    print("\n--- Another Example ---")
    # For an equilateral triangle a=b=c=s
    side_s = 10
    calculate_bm_mi_ratio(side_s, side_s, side_s)