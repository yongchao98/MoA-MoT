import math

def solve_lissajous_intersections():
    """
    Calculates the number of self-intersections for a Lissajous curve
    of the form x=cos(at), y=sin(bt).
    """
    # The curve is defined by (cos(9t), sin(5t))
    # This corresponds to a Lissajous curve with parameters a and b.
    a = 9
    b = 5

    # The formula for the number of self-intersections is (a-1)*(b-1)/2,
    # provided that a and b are coprime integers.
    # First, we check if a and b are coprime using the Greatest Common Divisor (GCD).
    gcd_val = math.gcd(a, b)

    if gcd_val == 1:
        # Calculate the intermediate values for the equation step-by-step
        a_minus_1 = a - 1
        b_minus_1 = b - 1
        numerator = a_minus_1 * b_minus_1
        result = numerator // 2

        # Print the final equation with all its numbers
        print(f"For the curve (cos({a}t), sin({b}t)), the number of self-intersections is found using the formula (a - 1) * (b - 1) / 2.")
        print("Calculation:")
        print(f"({a} - 1) * ({b} - 1) / 2 = {a_minus_1} * {b_minus_1} / 2 = {numerator} / 2 = {result}")
    else:
        print(f"The standard formula cannot be used because a={a} and b={b} are not coprime (GCD = {gcd_val}).")

solve_lissajous_intersections()
<<<16>>>