import math

def find_lissajous_intersections():
    """
    Calculates and explains the number of self-intersection points for the
    Lissajous curve defined by (cos(9t), sin(5t)).
    """
    # Define the parameters from the curve equations
    a = 9
    b = 5

    # The formula for the number of self-intersections of a Lissajous curve
    # (cos(at), sin(bt)) for coprime integers a and b is (a-1)(b-1)/2.
    
    # Check if a and b are coprime
    if math.gcd(a, b) == 1:
        # Calculate the components of the formula
        a_minus_1 = a - 1
        b_minus_1 = b - 1
        numerator = a_minus_1 * b_minus_1
        num_intersections = numerator // 2

        print("The number of self-intersection points for the curve (cos(at), sin(bt)) is calculated using the formula:")
        print("N = (a - 1) * (b - 1) / 2")
        print(f"\nFor the given curve, a = {a} and b = {b}.")
        print("The calculation is as follows:")
        print(f"N = ({a} - 1) * ({b} - 1) / 2")
        print(f"N = {a_minus_1} * {b_minus_1} / 2")
        print(f"N = {numerator} / 2")
        print(f"N = {num_intersections}")
    else:
        print(f"The parameters a={a} and b={b} are not coprime, so the standard formula does not apply.")

# Run the calculation
find_lissajous_intersections()