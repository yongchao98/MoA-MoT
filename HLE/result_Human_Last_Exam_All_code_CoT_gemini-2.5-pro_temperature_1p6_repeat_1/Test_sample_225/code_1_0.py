import math

def solve():
    """
    Calculates the number of self-intersection points for the curve (cos(9t), sin(5t)).
    """
    # The curve is defined by x = cos(at) and y = sin(bt)
    a = 9
    b = 5

    print("The given curve is a Lissajous curve of the form (cos(a*t), sin(b*t)).")
    print(f"For this curve, the parameters are a = {a} and b = {b}.")
    
    # Check if a and b are coprime and both odd
    are_coprime = math.gcd(a, b) == 1
    are_both_odd = a % 2 != 0 and b % 2 != 0

    if are_coprime and are_both_odd:
        print("\nSince a and b are coprime and both are odd, we can use the formula for self-intersections:")
        print("Number of self-intersections = (a - 1) * (b - 1) / 2")

        # Calculate the intermediate values
        a_minus_1 = a - 1
        b_minus_1 = b - 1
        numerator = a_minus_1 * b_minus_1
        
        # Calculate the final result
        result = numerator // 2

        print("\nSubstituting the values into the formula:")
        # The task requires printing each number in the final equation.
        # So we construct the string carefully.
        print(f"Number of self-intersections = ({a} - 1) * ({b} - 1) / 2")
        print(f"                             = {a_minus_1} * {b_minus_1} / 2")
        print(f"                             = {numerator} / 2")
        print(f"                             = {result}")
    else:
        print("The standard formula does not apply because the conditions are not met.")

solve()