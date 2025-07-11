import math

def solve_lissajous_intersections():
    """
    Calculates the number of self-intersection points for the Lissajous curve
    defined by x = cos(a*t) and y = sin(b*t).
    """
    # For the curve (cos(9t), sin(5t)), the parameters are a=9 and b=5.
    a = 9
    b = 5

    print(f"The parametric curve is (cos({a}t), sin({b}t)), which is a Lissajous curve.")
    print(f"The parameters are a = {a} and b = {b}.")

    # The formula for the number of self-intersections of a Lissajous curve is (a-1)(b-1)/2,
    # provided that a and b are coprime (their greatest common divisor is 1).

    # First, we check if a and b are coprime.
    common_divisor = math.gcd(a, b)
    if common_divisor == 1:
        print(f"The greatest common divisor of {a} and {b} is {common_divisor}, so they are coprime.")
        
        # Calculate the components of the formula.
        term1 = a - 1
        term2 = b - 1
        numerator = term1 * term2
        num_intersections = numerator // 2

        # Print the final calculation step-by-step.
        print("\nApplying the formula for the number of self-intersections: (a-1)*(b-1)/2")
        print(f"({a} - 1) * ({b} - 1) / 2 = {term1} * {term2} / 2 = {numerator} / 2 = {num_intersections}")

        print(f"\nTherefore, the number of self-intersection points is {num_intersections}.")
    else:
        print(f"The parameters a={a} and b={b} are not coprime (GCD={common_divisor}), so the standard formula does not apply.")

solve_lissajous_intersections()