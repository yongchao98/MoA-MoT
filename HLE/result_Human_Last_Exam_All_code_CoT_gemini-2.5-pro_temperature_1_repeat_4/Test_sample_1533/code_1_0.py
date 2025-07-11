def calculate_bm_mi_ratio(a, b, c):
    """
    Calculates the ratio BM/MI in a triangle ABC.

    Args:
        a: The length of side BC, opposite to vertex A.
        b: The length of side AC, opposite to vertex B.
        c: The length of side AB, opposite to vertex C.

    Returns:
        The numerical value of the ratio (a+c)/b, or an error string if the
        side lengths do not form a valid triangle.
    """
    # The triangle inequality theorem states that the sum of any two sides
    # of a triangle must be greater than the third side.
    if a + b <= c or a + c <= b or b + c <= a:
        return "Error: The given side lengths do not form a valid triangle."

    # The ratio BM/MI is given by the formula (a+c)/b.
    ratio = (a + c) / b
    return ratio

# --- Example Usage ---
# Let's consider a triangle with side lengths a=7, b=8, c=9.
# a is the side opposite to angle A (BC)
# b is the side opposite to angle B (AC)
# c is the side opposite to angle C (AB)
side_a = 7
side_b = 8
side_c = 9

# Calculate the ratio for the example triangle
result = calculate_bm_mi_ratio(side_a, side_b, side_c)

# Print the final result and the equation with the numbers used
if isinstance(result, str):
    print(result)
else:
    print(f"For a triangle with side lengths a = {side_a}, b = {side_b}, c = {side_c}:")
    print("The ratio BM / MI is given by the formula (a + c) / b.")
    print(f"BM / MI = ({side_a} + {side_c}) / {side_b} = {result}")
