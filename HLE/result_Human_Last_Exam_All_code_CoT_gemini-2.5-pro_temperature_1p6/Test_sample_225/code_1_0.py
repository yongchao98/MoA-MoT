import math

def solve():
    """
    This function calculates the number of self-intersection points
    for a Lissajous curve with parameters a=9 and b=5.
    """
    # Parameters from the curve equations x=cos(at), y=sin(bt)
    a = 9
    b = 5

    # Check if the conditions for the formula are met.
    # a and b should be coprime and both odd.
    if math.gcd(a, b) == 1 and a % 2 != 0 and b % 2 != 0:
        # The formula for the number of self-intersections for a Lissajous curve
        # where a and b are odd coprime integers.
        intersections = (a - 1) * (b - 1) // 2

        print(f"The curve is a Lissajous curve with parameters a = {a} and b = {b}.")
        print("The number of self-intersection points is found using the formula:")
        print("N = (a - 1) * (b - 1) / 2")
        print(f"N = ({a} - 1) * ({b} - 1) / 2")
        print(f"N = {a-1} * {b-1} / 2")
        print(f"N = {(a-1)*(b-1)} / 2")
        print(f"N = {intersections}")
    else:
        print("The standard formula does not apply for the given a and b.")

solve()