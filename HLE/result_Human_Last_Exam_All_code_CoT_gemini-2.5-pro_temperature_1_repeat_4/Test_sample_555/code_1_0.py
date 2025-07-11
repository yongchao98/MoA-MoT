import sympy

def solve_expected_rolls():
    """
    Calculates the minimal expected number of rolls of a 7-sided die
    to generate a uniform random digit (0-9).
    """
    # Define the symbols for the expected values E_s, where s is the number of states.
    # The states involved in the cycle are 1, 7, 49, 9, 63, 3, 21.
    E1, E3, E7, E9, E21, E49, E63 = sympy.symbols('E1 E3 E7 E9 E21 E49 E63')

    # Define the system of equations based on the optimal strategy.
    # If states < 10, E_s = 1 + E_{s*7}
    # If states >= 10, E_s = (s % 10 / s) * E_{s % 10}
    
    # Equation 1: Start with 1 state. Must roll.
    eq1 = sympy.Eq(E1, 1 + E7)
    
    # Equation 2: From 7 states. Must roll.
    eq2 = sympy.Eq(E7, 1 + E49)
    
    # Equation 3: From 49 states. Can decide. 9 states are leftovers.
    eq3 = sympy.Eq(E49, (49 % 10) / 49 * E9)
    
    # Equation 4: From 9 states. Must roll.
    eq4 = sympy.Eq(E9, 1 + E63)
    
    # Equation 5: From 63 states. Can decide. 3 states are leftovers.
    eq5 = sympy.Eq(E63, (63 % 10) / 63 * E3)

    # Equation 6: From 3 states. Must roll.
    eq6 = sympy.Eq(E3, 1 + E21)

    # Equation 7: From 21 states. Can decide. 1 state is leftover, returning to the start.
    eq7 = sympy.Eq(E21, (21 % 10) / 21 * E1)

    # Solve the system of linear equations for all E_s variables.
    solution = sympy.solve([eq1, eq2, eq3, eq4, eq5, eq6, eq7], (E1, E3, E7, E9, E21, E49, E63))

    # The result for E1 is the answer we are looking for.
    # sympy automatically provides it as a simplified fraction.
    e1_value = solution[E1]
    
    # Extract numerator and denominator for the final output.
    numerator = e1_value.p
    denominator = e1_value.q

    print(f"The minimal expected value of rolls, E, is found by solving a system of linear equations.")
    print(f"The final equation for the expected value is:")
    print(f"E = {numerator} / {denominator}")
    
if __name__ == '__main__':
    solve_expected_rolls()