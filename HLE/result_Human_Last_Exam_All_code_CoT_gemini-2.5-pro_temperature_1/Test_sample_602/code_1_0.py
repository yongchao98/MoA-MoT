import math

def solve_l_n():
    """
    This function calculates the exact value of l(n) as n -> infinity.

    The derived expression for l(n) can be written in the form:
    l(n) = C1 + (C2 + (C3*n + C4) * sqrt(C5*n^2 + C6*n + C7)) / n**C8
    """

    # Coefficients of the derived formula for l(n)
    # l(n) = 2 + (2 + (-4*n + 2) * sqrt(1*n^2 - 1*n + 1)) / n**2
    C1 = 2
    C2 = 2
    C3 = -4
    C4 = 2
    C5 = 1
    C6 = -1
    C7 = 1
    C8 = 2

    print("The derived expression for l(n) has the form:")
    print("l(n) = C1 + (C2 + (C3*n + C4) * sqrt(C5*n^2 + C6*n + C7)) / n**C8\n")

    print("The integer coefficients are:")
    print(f"C1 = {C1}")
    print(f"C2 = {C2}")
    print(f"C3 = {C3}")
    print(f"C4 = {C4}")
    print(f"C5 = {C5}")
    print(f"C6 = {C6}")
    print(f"C7 = {C7}")
    print(f"C8 = {C8}\n")

    # The exact value is the limit of l(n) as n -> infinity.
    # The limit can be calculated from the coefficients as: C1 + C3 * sqrt(C5)
    limit_value = C1 + C3 * math.sqrt(C5)

    print("The exact value of l(n) is its limit as n approaches infinity.")
    print(f"Limit = C1 + C3 * sqrt(C5) = {C1} + {C3} * sqrt({C5}) = {limit_value}")

if __name__ == '__main__':
    solve_l_n()