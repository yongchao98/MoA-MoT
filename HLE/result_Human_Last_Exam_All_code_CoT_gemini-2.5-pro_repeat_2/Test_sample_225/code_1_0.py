import math

def solve():
    """
    Calculates the number of self-intersection points for the curve
    (cos(9t), sin(5t)).
    """
    # The curve is a Lissajous curve of the form (cos(at), sin(bt)).
    a = 9
    b = 5

    # Check if a and b are coprime
    if math.gcd(a, b) != 1:
        print(f"The numbers a={a} and b={b} are not coprime.")
        print("The standard formula may not apply directly.")
        return

    # The formula for the number of self-intersections for a Lissajous curve
    # of the form (cos(at), sin(bt)) with a and b coprime is (a-1)(b-1)/2.
    a_minus_1 = a - 1
    b_minus_1 = b - 1
    numerator = a_minus_1 * b_minus_1
    num_intersections = numerator // 2

    print("The curve is a Lissajous curve of the form (cos(at), sin(bt)).")
    print(f"Here, a = {a} and b = {b}.")
    print("The formula for the number of self-intersection points is (a - 1) * (b - 1) / 2.")
    print("\nCalculation:")
    print(f"Number of intersections = ({a} - 1) * ({b} - 1) / 2")
    print(f"                         = {a_minus_1} * {b_minus_1} / 2")
    print(f"                         = {numerator} / 2")
    print(f"                         = {num_intersections}")

solve()
