import math

def solve_self_intersections():
    """
    Calculates the number of self-intersection points for the Lissajous curve
    defined by x(t) = cos(at) and y(t) = sin(bt).
    """
    # Parameters from the curve equations x(t) = cos(9t), y(t) = sin(5t)
    a = 9
    b = 5

    # The formula for the number of self-intersections for a Lissajous curve
    # of this type (where a and b are odd) is (a * b - gcd(a, b)) / 2.
    
    # First, calculate the greatest common divisor (GCD) of a and b.
    gcd_val = math.gcd(a, b)

    # Check if the conditions for the formula are met (optional, but good practice)
    # This specific formula works when a and b are odd integers.
    if a % 2 == 1 and b % 2 == 1:
        # Apply the formula to find the number of self-intersections.
        num_intersections = (a * b - gcd_val) // 2

        # Print the equation with the specific numbers, as requested.
        print(f"The number of self-intersections is calculated using the formula (a * b - gcd(a, b)) / 2.")
        print(f"For a = {a} and b = {b}:")
        print(f"({a} * {b} - {gcd_val}) / 2 = {num_intersections}")
    else:
        print("The formula used applies to cases where a and b are odd integers.")

solve_self_intersections()