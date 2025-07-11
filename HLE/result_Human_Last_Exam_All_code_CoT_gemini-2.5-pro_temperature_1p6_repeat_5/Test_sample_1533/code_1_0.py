def calculate_bm_mi_ratio(a, b, c):
    """
    Calculates the ratio BM/MI for a triangle ABC with side lengths a, b, c.

    In a triangle ABC:
    - a is the length of side BC (opposite to angle A)
    - b is the length of side CA (opposite to angle B)
    - c is the length of side AB (opposite to angle C)
    - I is the incenter.
    - M is the intersection of the angle bisector BI with the circumcircle.

    The ratio BM/MI is given by the formula (a + c) / b.

    Args:
    a (float): Length of side BC.
    b (float): Length of side CA.
    c (float): Length of side AB.

    Returns:
    None. Prints the result.
    """
    if a + b <= c or a + c <= b or b + c <= a:
        print("The given side lengths do not form a valid triangle.")
        return

    # The ratio is derived to be (a + c) / b
    ratio_numerator = a + c
    ratio_denominator = b
    ratio_value = ratio_numerator / ratio_denominator

    # Print the formula and the calculated values
    print("The ratio BM/MI is given by the formula (a + c) / b.")
    print(f"For a triangle with side lengths a = {a}, b = {b}, c = {c}:")
    print(f"The numerator of the ratio is a + c = {a} + {c} = {ratio_numerator}")
    print(f"The denominator of the ratio is b = {b}")
    print(f"Therefore, the final ratio BM/MI = {ratio_numerator} / {b} = {ratio_value}")


# Example usage with a 3-4-5 right triangle where b=4
# Sides are a=BC=3, b=CA=4, c=AB=5. Angle C is 90 degrees.
calculate_bm_mi_ratio(a=3, b=4, c=5)