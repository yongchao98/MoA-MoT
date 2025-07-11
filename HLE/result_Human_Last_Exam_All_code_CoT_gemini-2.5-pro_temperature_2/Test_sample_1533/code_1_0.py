def calculate_bm_mi_ratio(a, b, c):
    """
    Calculates the ratio BM/MI in a triangle ABC.

    Args:
        a (float): The length of side BC, opposite to vertex A.
        b (float): The length of side AC, opposite to vertex B.
        c (float): The length of side AB, opposite to vertex C.
    """
    # Check if the given side lengths can form a triangle (Triangle Inequality Theorem)
    if (a + b <= c) or (a + c <= b) or (b + c <= a):
        print(f"The side lengths {a}, {b}, and {c} cannot form a triangle.")
        return

    # The derived formula for the ratio BM/MI is (a + c) / b
    ratio = (a + c) / b

    # Print the explanation and the result
    print("In a triangle ABC with side lengths a, b, and c:")
    print("The angle bisector of ∠ABC (BI) intersects the circumcircle at M.")
    print("The ratio BM/MI is given by the formula: (a + c) / b")
    print("\nFor the given triangle with sides a={}, b={}, c={}:".format(a, b, c))
    print("The final equation is: ({} + {}) / {}".format(a, c, b))
    print("Resulting ratio BM/MI = {}".format(ratio))

# Example usage with a common 3-4-5 right triangle
# Here, a=3, b=4, c=5
print("Example 1: a=3, b=4, c=5")
calculate_bm_mi_ratio(3, 4, 5)

print("-" * 30)

# Example usage with an equilateral triangle
# Here, a=6, b=6, c=6
print("Example 2: a=6, b=6, c=6")
calculate_bm_mi_ratio(6, 6, 6)

print("-" * 30)

# Example usage with an isosceles triangle
# Here, a=7, b=10, c=7
print("Example 3: a=7, b=10, c=7")
calculate_bm_mi_ratio(7, 10, 7)