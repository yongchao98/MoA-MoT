from fractions import Fraction

def solve_for_expected_rolls():
    """
    This function solves the system of linear equations for the expected number of rolls.

    Let E_C be the expected number of additional rolls required, given C equiprobable outcomes.
    The system of equations is:
    E_1 = 1 + E_7
    E_7 = 1 + (9/49) * E_9
    E_9 = 1 + (1/21) * E_3
    E_3 = 1 + (1/21) * E_1
    
    We solve this system to find E_1.
    """
    
    # E_3 in terms of E_1
    # E_3 = 1 + (1/21) * E_1
    c3 = Fraction(1)
    d3 = Fraction(1, 21)
    
    # E_9 in terms of E_3, then substitute to get it in terms of E_1
    # E_9 = 1 + (1/21) * E_3 = 1 + (1/21) * (c3 + d3*E_1)
    c9_from_e3 = Fraction(1)
    d9_from_e3 = Fraction(1, 21)
    c9 = c9_from_e3 + d9_from_e3 * c3
    d9 = d9_from_e3 * d3
    
    # E_7 in terms of E_9, then substitute to get it in terms of E_1
    # E_7 = 1 + (9/49) * E_9
    c7_from_e9 = Fraction(1)
    d7_from_e9 = Fraction(9, 49)
    c7 = c7_from_e9 + d7_from_e9 * c9
    d7 = d7_from_e9 * d9
    
    # E_1 in terms of E_7, then substitute to get an equation for E_1
    # E_1 = 1 + E_7
    c1_from_e7 = Fraction(1)
    d1_from_e7 = Fraction(1)
    c1_final = c1_from_e7 + d1_from_e7 * c7
    d1_final = d1_from_e7 * d7

    # We now have an equation of the form: E_1 = c1_final + d1_final * E_1
    # E_1 * (1 - d1_final) = c1_final
    # E_1 = c1_final / (1 - d1_final)
    
    e1_numerator = c1_final
    e1_denominator = 1 - d1_final
    
    final_result = e1_numerator / e1_denominator
    
    print("The minimal expected value E is found by solving a system of linear equations.")
    print("The final equation for the expected number of rolls, E, is of the form: E = C + D * E")
    print(f"In our case, we derived E = {c1_final.numerator}/{c1_final.denominator} + ({d1_final.numerator}/{d1_final.denominator}) * E")
    print("Solving for E gives E = C / (1 - D):")
    
    intermediate_calc_1 = f"E = ({c1_final.numerator}/{c1_final.denominator}) / (1 - {d1_final.numerator}/{d1_final.denominator})"
    intermediate_calc_2 = f"E = ({e1_numerator.numerator}/{e1_numerator.denominator}) / ({e1_denominator.numerator}/{e1_denominator.denominator})"
    
    # Let's perform the calculation in a more familiar form for display
    # E = c1 / (1-d1) = (752/343) / (2400/2401) = (752/343) * (2401/2400)
    display_calc = f"E = ({e1_numerator.numerator}/{e1_numerator.denominator}) * ({e1_denominator.denominator}/{e1_denominator.numerator})"
    unsimplified_fraction = e1_numerator * e1_denominator.denominator
    unsimplified_denominator = e1_denominator.numerator
    
    
    print(intermediate_calc_1)
    print(intermediate_calc_2)
    print(display_calc)
    print(f"E = {unsimplified_fraction}/{unsimplified_denominator}")
    print(f"The simplified result is: {final_result.numerator}/{final_result.denominator}")


solve_for_expected_rolls()