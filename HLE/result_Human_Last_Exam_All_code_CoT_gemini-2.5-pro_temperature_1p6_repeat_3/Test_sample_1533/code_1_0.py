def solve_ratio(a, b, c):
    """
    Calculates the ratio BM/MI for a triangle ABC with side lengths a, b, c.

    In triangle ABC:
    - a is the length of side BC (opposite angle A)
    - b is the length of side AC (opposite angle B)
    - c is the length of side AB (opposite angle C)
    - I is the incenter.
    - M is the intersection of angle bisector BI with the circumcircle.

    The ratio BM/MI is given by the formula (a + c) / b.
    """
    # Check if the side lengths can form a triangle
    if a + b <= c or a + c <= b or b + c <= a:
        print("The given side lengths do not form a valid triangle.")
        return

    # The derived ratio is (a + c) / b
    result = (a + c) / b

    print("The problem asks for the ratio BM/MI.")
    print("The derived expression in terms of side lengths a, b, and c is: (a + c) / b")
    print("\nFor the triangle with side lengths:")
    print(f"a = {a}")
    print(f"b = {b}")
    print(f"c = {c}")
    
    print("\nThe calculation is:")
    # Per the instruction, outputting each number in the final equation
    print(f"Ratio = ({a} + {c}) / {b}")
    
    print("\nFinal Result:")
    print(f"The ratio BM/MI is {result}")

# Example usage with a common scalene triangle (sides 7, 8, 9)
# Let a=7 (opposite A), b=8 (opposite B), c=9 (opposite C)
solve_ratio(a=7, b=8, c=9)
