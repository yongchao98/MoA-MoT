def calculate_bm_mi_ratio(a, b, c):
    """
    Calculates the ratio BM/MI in a triangle ABC.

    The problem defines:
    - a, b, c as the side lengths of triangle ABC.
    - I as the incenter.
    - M as the intersection of the angle bisector BI and the circumcircle.

    The derived formula for the ratio BM/MI is (a + c) / b.

    Args:
      a: The length of side BC (opposite vertex A).
      b: The length of side AC (opposite vertex B).
      c: The length of side AB (opposite vertex C).
    """
    # First, verify if the given side lengths can form a valid triangle.
    if not (a + b > c and a + c > b and b + c > a):
        print(f"Error: The side lengths a={a}, b={b}, c={c} cannot form a triangle.")
        return

    # The ratio is (a+c) / b based on the geometric derivation.
    numerator = a + c
    ratio = numerator / b

    print("The ratio BM/MI is expressed by the formula: (a + c) / b")
    print(f"For a triangle with side lengths a = {a}, b = {b}, c = {c}:")
    print(f"BM / MI = ({a} + {c}) / {b} = {numerator} / {b} = {ratio}")

# Example usage with a sample triangle with sides a=7, b=8, c=9
calculate_bm_mi_ratio(7, 8, 9)