def solve_triangle_ratio(a, b, c):
    """
    This function calculates the ratio BM/MI based on the derived formula.

    In triangle ABC, with side lengths a, b, and c:
    - a is the length of side BC
    - b is the length of side CA
    - c is the length of side AB
    - I is the incenter.
    - M is the point where the angle bisector of angle B intersects the circumcircle.

    The ratio BM/MI is given by the formula (a + c) / b.
    This function takes the side lengths as input, checks if they can form a valid
    triangle, and prints the result of the formula.
    """
    
    # The triangle inequality theorem states that the sum of the lengths of any two sides
    # of a triangle must be greater than the length of the third side.
    if not (a + b > c and a + c > b and b + c > a):
        print(f"Error: The side lengths a={a}, b={b}, c={c} cannot form a valid triangle.")
        return
        
    # Calculate the ratio using the derived formula
    ratio = (a + c) / b
    
    print(f"For a triangle with side lengths a = {a}, b = {b}, c = {c}:")
    print("The ratio BM/MI is given by the formula (a + c) / b.")
    print(f"BM/MI = ({a} + {c}) / {b} = {ratio}")

# You can test this with any valid triangle side lengths.
# Let's use the common 3-4-5 right triangle as an example.
# Let a = 3 (BC), b = 4 (CA), and c = 5 (AB).
side_a = 3
side_b = 4
side_c = 5

solve_triangle_ratio(side_a, side_b, side_c)
