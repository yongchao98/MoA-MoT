def solve_lissajous_intersections():
    """
    Calculates the number of self-intersection points for the curve
    (cos(9t), sin(5t)).
    """
    # The parameters from the curve equations x=cos(at), y=sin(bt)
    a = 9
    b = 5

    # The formula for the number of self-intersections of a Lissajous curve
    # with coprime parameters a and b is N = (a - 1) * (b - 1) / 2.

    # Intermediate calculation steps
    a_minus_1 = a - 1
    b_minus_1 = b - 1
    numerator = a_minus_1 * b_minus_1
    num_intersections = numerator // 2

    # Print the explanation and the step-by-step calculation
    print("The curve is a Lissajous curve with parameters a = 9 and b = 5.")
    print("The number of self-intersections (N) is given by the formula: N = (a - 1) * (b - 1) / 2")
    print("\nCalculation steps:")
    print(f"N = ({a} - 1) * ({b} - 1) / 2")
    print(f"N = {a_minus_1} * {b_minus_1} / 2")
    print(f"N = {numerator} / 2")
    print(f"N = {num_intersections}")
    print(f"\nThe total number of self-intersection points is {num_intersections}.")

solve_lissajous_intersections()