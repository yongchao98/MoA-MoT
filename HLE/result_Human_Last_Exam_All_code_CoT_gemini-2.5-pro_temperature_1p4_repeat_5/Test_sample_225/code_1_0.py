def solve_self_intersections():
    """
    Calculates the number of self-intersection points for the curve
    (cos(9t), sin(5t)).
    """
    a = 9
    b = 5

    print("The curve is a Lissajous curve defined by x = cos(a*t) and y = sin(b*t).")
    print(f"For the given curve, the parameters are a = {a} and b = {b}.")
    print("\nTo find the number of self-intersections, we first check if 'a' and 'b' are coprime.")
    
    # Simple GCD check for explanation
    import math
    if math.gcd(a, b) == 1:
        print(f"The greatest common divisor of {a} and {b} is 1, so they are coprime.")
    else:
        print(f"The parameters a={a} and b={b} are not coprime. The formula might differ.")
        return

    print("\nBoth a=9 and b=5 are odd numbers.")
    print("The formula for the number of self-intersection points (N) for a Lissajous curve with coprime odd 'a' and 'b' is:")
    print("N = (a - 1) * (b - 1) / 2")

    print("\nPlugging in the values a = 9 and b = 5:")
    
    # Step-by-step calculation
    a_minus_1 = a - 1
    b_minus_1 = b - 1
    numerator = a_minus_1 * b_minus_1
    result = numerator / 2

    print(f"N = ({a} - 1) * ({b} - 1) / 2")
    print(f"N = {a_minus_1} * {b_minus_1} / 2")
    print(f"N = {numerator} / 2")
    print(f"N = {int(result)}")
    
    print(f"\nThe number of self-intersection points for the curve is {int(result)}.")

solve_self_intersections()