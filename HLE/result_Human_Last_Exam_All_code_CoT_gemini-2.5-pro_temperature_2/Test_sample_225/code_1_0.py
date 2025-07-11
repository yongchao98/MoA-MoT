import math

def solve():
    """
    Calculates the number of self-intersection points for the curve
    (cos(9t), sin(5t)).
    """
    # The curve is given by x(t) = cos(nt) and y(t) = sin(kt)
    # with n = 9 and k = 5.
    n = 9
    k = 5

    # This is a Lissajous curve. To use a standard formula, we write it as:
    # x(t) = cos(nt)
    # y(t) = cos(kt + psi), where sin(5t) = cos(5t - pi/2).

    # The number of self-intersections for coprime n and k is given by a formula.
    # First, confirm they are coprime.
    if math.gcd(n, k) != 1:
        print(f"Error: n={n} and k={k} are not coprime.")
        return

    # For n, k coprime and both odd, the formula is:
    # N = (n*k - n - k + 1) / 2
    
    # Calculate the number of intersections
    numerator_val = n * k - n - k + 1
    num_intersections = numerator_val / 2

    print("To find the number of self-intersection points for the curve (cos(9t), sin(5t)), we use a standard formula for Lissajous curves.")
    print("The parameters are n=9 and k=5.")
    print("Since n and k are coprime and both odd, the number of self-intersections is given by:")
    print("N = (n*k - n - k + 1) / 2")
    print(f"N = ({n}*{k} - {n} - {k} + 1) / 2")
    print(f"N = ({n*k} - {n+k} + 1) / 2")
    print(f"N = ({numerator_val}) / 2")
    print(f"N = {int(num_intersections)}")

solve()