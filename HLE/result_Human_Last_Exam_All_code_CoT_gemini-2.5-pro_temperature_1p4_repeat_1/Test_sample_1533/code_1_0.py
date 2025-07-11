import sympy

def solve_triangle_ratio_incenter():
    """
    This function provides the solution to the specified geometry problem by deriving
    the ratio BM/MI in terms of the side lengths a, b, and c of triangle ABC.

    The derivation relies on Ptolemy's Theorem for the cyclic quadrilateral ABCM
    and the properties of the point M, which is the intersection of an angle
    bisector with the circumcircle. The key properties are AM = CM and AM = MI.

    The final derived ratio is (a + c) / b.
    """

    # Define symbolic variables for the side lengths of the triangle ABC
    a, b, c = sympy.symbols('a b c')

    # According to the derivation, the numerator of the ratio is the sum of
    # the lengths of the sides adjacent to angle B.
    numerator_part_1 = a  # Side BC
    numerator_part_2 = c  # Side AB

    # The denominator is the length of the side opposite to angle B.
    denominator_part = b  # Side AC

    # The final equation for the ratio BM/MI is (a + c) / b.
    # The following print statement displays this result clearly,
    # showing each component of the final expression.
    print("The final expression for the ratio BM/MI is:")
    print(f"({numerator_part_1} + {numerator_part_2}) / {denominator_part}")

# Execute the function to display the solution
solve_triangle_ratio_incenter()