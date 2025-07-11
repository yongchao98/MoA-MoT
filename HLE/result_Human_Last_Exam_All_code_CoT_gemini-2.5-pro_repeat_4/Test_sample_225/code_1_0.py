import math

def solve():
    """
    Calculates the number of self-intersection points for the curve
    (cos(9t), sin(5t)).
    """
    # The curve is a Lissajous curve of the form x(t) = cos(at), y(t) = sin(bt).
    # The parameters for the given curve are:
    a = 9
    b = 5

    # The number of self-intersection points for such a curve, where 'a' and 'b'
    # are coprime odd integers, is given by the formula (a - 1) * (b - 1).
    # We can verify that a=9 and b=5 are indeed coprime odd integers.
    # math.gcd(9, 5) is 1.
    # 9 % 2 is 1 and 5 % 2 is 1.

    # Calculate the number of self-intersection points using the formula.
    num_intersections = (a - 1) * (b - 1)
    
    # Print the explanation and the calculation steps.
    print(f"The curve is (cos({a}t), sin({b}t)).")
    print(f"The number of self-intersection points is given by the formula (a-1) * (b-1).")
    print(f"Substituting the values a={a} and b={b}:")
    print(f"({a} - 1) * ({b} - 1) = {a - 1} * {b - 1} = {num_intersections}")
    print(f"The number of self-intersection points is {num_intersections}.")

solve()