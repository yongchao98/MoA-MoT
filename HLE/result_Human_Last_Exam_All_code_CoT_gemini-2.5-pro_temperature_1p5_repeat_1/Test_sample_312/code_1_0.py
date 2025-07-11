from fractions import Fraction

def solve_exponent():
    """
    Calculates the smallest possible value for the exponent c based on the
    dimension delta of the Frostman measure.
    """
    delta_num = 8
    delta_den = 5
    delta = Fraction(delta_num, delta_den)

    # The formula for the exponent c is (1 - delta) / 2
    # We will show the calculation step-by-step.

    # Step 1: Start with the formula for c
    print(f"The smallest possible value for c is given by the formula:")
    print(f"c = (1 - delta) / 2")

    # Step 2: Substitute the value of delta
    print(f"\nGiven delta = {delta_num}/{delta_den}, we substitute this into the formula:")
    print(f"c = (1 - {delta_num}/{delta_den}) / 2")

    # Step 3: Calculate the term in the parenthesis
    numerator_in_parenthesis = delta_den - delta_num
    print(f"\nFirst, calculate the expression in the parenthesis:")
    print(f"1 - {delta_num}/{delta_den} = {delta_den}/{delta_den} - {delta_num}/{delta_den} = {numerator_in_parenthesis}/{delta_den}")

    # Step 4: Substitute back and perform the final division
    c = Fraction(numerator_in_parenthesis, delta_den) / 2
    print(f"\nNow substitute this back into the equation for c:")
    print(f"c = ({numerator_in_parenthesis}/{delta_den}) / 2")
    print(f"c = {c.numerator}/{c.denominator}")

solve_exponent()