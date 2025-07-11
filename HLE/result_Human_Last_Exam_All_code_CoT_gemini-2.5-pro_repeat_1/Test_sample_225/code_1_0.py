import math

def find_self_intersections():
    """
    Calculates the number of self-intersection points for the curve
    defined by x = cos(at) and y = sin(bt).
    """
    # The parameters from the curve equations x = cos(9t), y = sin(5t)
    a = 9
    b = 5

    print(f"The curve is a Lissajous curve with parameters a = {a} and b = {b}.")
    
    # The formula for the number of self-intersections is (a-1)(b-1)/2,
    # given that a and b are odd and coprime.
    
    # Verify conditions
    is_coprime = math.gcd(a, b) == 1
    a_is_odd = a % 2 != 0
    b_is_odd = b % 2 != 0

    if is_coprime and a_is_odd and b_is_odd:
        print("Conditions met: 'a' and 'b' are odd and coprime.")
        
        # Calculate the number of self-intersections
        numerator = (a - 1) * (b - 1)
        num_intersections = numerator // 2
        
        print("\nCalculating the number of self-intersections using the formula N = (a - 1) * (b - 1) / 2:")
        print(f"N = ({a} - 1) * ({b} - 1) / 2")
        print(f"N = {a - 1} * {b - 1} / 2")
        print(f"N = {numerator} / 2")
        print(f"N = {num_intersections}")
        
        print(f"\nThe total number of self-intersection points is {num_intersections}.")
    else:
        print("The conditions for the simple formula are not met.")

find_self_intersections()