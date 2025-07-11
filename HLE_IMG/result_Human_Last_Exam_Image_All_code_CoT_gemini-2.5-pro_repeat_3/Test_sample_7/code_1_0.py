def get_plateau_formula():
    """
    This function determines and prints the formula for the second voltage plateau
    of a graphite anode during lithium intercalation.
    """

    # The formula relates voltage (V) to chemical potential (μ) and elementary charge (e).
    # For the second plateau, the relevant chemical potential is μ_2, which corresponds
    # to the formation of the Stage 2 structure.

    # Define the numerical components of the formula
    coefficient = -1
    stage_index = 2

    # Define the symbols
    potential_symbol = "μ"
    charge_symbol = "e"
    voltage_symbol = "V"

    # Print the formula, explicitly showing the numbers involved.
    print("The simple formula for the second voltage plateau is:")
    print(f"{voltage_symbol} = ({coefficient} * {potential_symbol}_{stage_index}) / {charge_symbol}")

    # Print the formula in its conventional, simplified form for clarity.
    print("\nOr, more simply:")
    print(f"{voltage_symbol} = -{potential_symbol}_{stage_index} / {charge_symbol}")

# Execute the function
get_plateau_formula()