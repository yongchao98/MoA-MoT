def solve_curve_intersections():
    """
    This function calculates the number of self-intersection points for a Lissajous curve
    defined by the parametric equations x(t) = cos(at) and y(t) = sin(bt).
    """
    # The curve is defined by (cos(9t), sin(5t)).
    # The parameters are a = 9 and b = 5.
    a = 9
    b = 5

    # The formula for the number of self-intersections for a Lissajous curve
    # with coprime integer parameters a and b is N = (a - 1) * (b - 1) / 2.
    # First, we check if a and b are coprime.
    def gcd(p, q):
        while q:
            p, q = q, p % q
        return p

    if gcd(a, b) != 1:
        print(f"The formula may not apply as a={a} and b={b} are not coprime.")
        return

    # Calculate the number of self-intersections.
    # The formula is N = (a - 1) * (b - 1) / 2
    part1 = a - 1
    part2 = b - 1
    numerator = part1 * part2
    num_intersections = numerator / 2

    # Print the explanation and the step-by-step calculation of the final equation.
    print(f"For the curve (cos({a}t), sin({b}t)), the number of self-intersections is calculated using the formula:")
    print("N = (a - 1) * (b - 1) / 2")
    print("\nSubstituting the values a = {} and b = {}:".format(a, b))
    print("N = ({} - 1) * ({} - 1) / 2".format(a, b))
    print("N = {} * {} / 2".format(part1, part2))
    print("N = {} / 2".format(numerator))
    print("N = {}".format(int(num_intersections)))

solve_curve_intersections()