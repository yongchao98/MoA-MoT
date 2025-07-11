import math

def solve_lissajous_intersections():
    """
    Calculates the number of self-intersection points for a Lissajous curve
    of the form (cos(at), sin(bt)) where a and b are odd coprime integers.
    """
    # Parameters from the curve equations x=cos(9t), y=sin(5t)
    a = 9
    b = 5

    # Check if a and b are coprime
    if math.gcd(a, b) != 1:
        print("The numbers a and b must be coprime.")
        return

    # Check if a and b are odd
    if a % 2 == 0 or b % 2 == 0:
        print("The numbers a and b must be odd for this specific formula.")
        return
    
    print("The curve is a Lissajous curve defined by (cos(at), sin(bt)).")
    print(f"For this curve, a = {a} and b = {b}.")
    print("Since a and b are coprime odd integers, the number of self-intersections is given by the formula (a - 1)(b - 1) / 2.")
    
    # Calculate the number of self-intersections
    intersections = (a - 1) * (b - 1) // 2
    
    # Print the full equation and the result
    print("\nCalculation:")
    print(f"({a} - 1) * ({b} - 1) / 2 = ({a-1}) * ({b-1}) / 2 = { (a-1)*(b-1) } / 2 = {intersections}")
    print(f"\nThe number of self-intersection points is {intersections}.")

solve_lissajous_intersections()
