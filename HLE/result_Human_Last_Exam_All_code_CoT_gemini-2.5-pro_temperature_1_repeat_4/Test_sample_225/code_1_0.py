import math

def solve_intersections():
    """
    Calculates the number of self-intersection points for the curve
    (cos(9t), sin(5t)).
    """
    # The curve is a Lissajous curve of the form x(t) = cos(nt), y(t) = sin(mt).
    # The parameters are n and m.
    n = 9
    m = 5

    # The number of self-intersection points for such a curve,
    # where n and m are coprime, is given by the formula:
    # N = (n - 1) * (m - 1) / 2
    
    # First, confirm that n and m are coprime.
    if math.gcd(n, m) != 1:
        print(f"Error: n={n} and m={m} are not coprime.")
        return

    # Apply the formula and show the steps of the calculation.
    n_minus_1 = n - 1
    m_minus_1 = m - 1
    numerator = n_minus_1 * m_minus_1
    num_intersections = numerator // 2

    print(f"Finding the number of self-intersections for the curve (cos({n}t), sin({m}t)).")
    print("The formula for the number of self-intersections is (n - 1) * (m - 1) / 2.")
    print("\nCalculation steps:")
    print(f"({n} - 1) * ({m} - 1) / 2")
    print(f"= {n_minus_1} * {m_minus_1} / 2")
    print(f"= {numerator} / 2")
    print(f"= {num_intersections}")
    
    print(f"\nThe total number of self-intersection points is {num_intersections}.")

solve_intersections()