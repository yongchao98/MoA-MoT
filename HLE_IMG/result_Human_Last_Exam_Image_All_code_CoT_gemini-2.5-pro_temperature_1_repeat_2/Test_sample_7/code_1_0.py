def solve_potential_formula():
    """
    This function determines and prints the formula for the second voltage plateau
    of a graphite anode during lithium intercalation.

    The voltage (V) of a plateau is determined by the chemical potential (μ)
    of the stage being formed, according to the relation V = -μ / e, where e
    is the elementary charge.

    The second plateau corresponds to the formation of Stage 2. The chemical
    potential for this process is denoted as μ_2.
    """

    # Define the components of the formula
    potential_symbol = "μ"
    stage_number = 2
    charge_symbol = "e"

    # Construct and print the formula for the second plateau
    # The equation includes the number 2 as the stage identifier.
    print(f"The formula for the second plateau is: V = -{potential_symbol}_{stage_number} / {charge_symbol}")

solve_potential_formula()