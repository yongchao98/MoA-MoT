import sympy

def solve_equation():
    """
    This function formulates and prints the simple formula for the second voltage plateau
    of a graphite anode based on the provided physical description.
    """

    # Define the numbers that will appear in the equation.
    # This fulfills the requirement to output each number in the final equation.
    stage_2_index = 2
    stage_3_index = 3
    denominator_value = 2

    # Define the symbols for the chemical potential (μ) and elementary charge (e).
    mu_symbol = "μ"
    e_symbol = "e"

    # Construct the formula string using an f-string.
    # The formula represents the average of the chemical potentials of the two
    # coexisting phases (Stage 2 and Stage 3) during the second plateau.
    formula = f"-({mu_symbol}_{stage_2_index} + {mu_symbol}_{stage_3_index}) / ({denominator_value}*e)"
    
    print("A simple formula that best approximates the second plateau is:")
    print(formula)

solve_equation()