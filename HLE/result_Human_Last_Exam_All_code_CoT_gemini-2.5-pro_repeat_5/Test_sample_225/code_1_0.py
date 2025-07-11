import math

def find_lissajous_intersections():
    """
    Calculates and explains the number of self-intersection points for the
    Lissajous curve (cos(9t), sin(5t)).
    """
    n = 9
    m = 5

    print(f"The problem is to find the number of self-intersection points for the curve (cos({n}t), sin({m}t)).")
    print("This is a Lissajous curve.")
    print("\nFor a Lissajous curve with parameters n and m that are coprime (gcd(n,m)=1),")
    print("the number of self-intersections is given by the formula:")
    print("N = (n - 1) * (m - 1) / 2\n")

    gcd_nm = math.gcd(n, m)
    print(f"First, let's confirm that n={n} and m={m} are coprime:")
    print(f"The greatest common divisor of {n} and {m} is {gcd_nm}.")
    if gcd_nm != 1:
        print("The parameters are not coprime. A different formula is needed.")
        return

    print("Since they are coprime, we can proceed with the formula.\n")

    # Step-by-step calculation
    n_minus_1 = n - 1
    m_minus_1 = m - 1
    numerator = n_minus_1 * m_minus_1
    result = numerator // 2

    print("Now we plug the numbers into the formula:")
    print(f"N = ({n} - 1) * ({m} - 1) / 2")
    print(f"N = {n_minus_1} * {m_minus_1} / 2")
    print(f"N = {numerator} / 2")
    print(f"N = {result}")

    print(f"\nThe number of self-intersection points is {result}.")

# Run the function to display the solution
find_lissajous_intersections()