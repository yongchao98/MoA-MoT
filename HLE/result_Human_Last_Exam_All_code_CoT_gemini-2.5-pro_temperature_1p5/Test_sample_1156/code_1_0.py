import math

def final_answer():
    """
    This function calculates and prints the normalised density of the invariant measure for the given map.
    """
    # The final equation for the normalised density rho(x) has been determined to be
    # rho(x) = (1/3) * (1 + 1/sqrt(x))

    # Define the numbers in the final equation
    coeff_numerator = 1
    coeff_denominator = 3
    
    term_1_constant = 1
    
    term_2_numerator = 1
    term_2_denominator = "sqrt(x)" # a symbolic representation

    print("The normalised density rho(x) for the map T(x) = 1/sqrt(x) mod 1 is:")
    print(f"rho(x) = ({coeff_numerator}/{coeff_denominator}) * ({term_1_constant} + {term_2_numerator}/{term_2_denominator})")

final_answer()