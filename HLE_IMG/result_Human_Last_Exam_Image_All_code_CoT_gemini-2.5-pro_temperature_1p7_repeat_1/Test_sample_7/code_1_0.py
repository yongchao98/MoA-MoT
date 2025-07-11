def print_voltage_formula():
    """
    This function prints the derived formula for the second voltage plateau
    of a graphite anode during lithium intercalation.
    """
    # Symbols used in the formula
    voltage_symbol = "V"
    mu_symbol = "μ"
    charge_symbol = "e"

    # Numeric indices and values in the formula
    stage_index_2 = 2
    stage_index_3 = 3
    denominator_coefficient = 2

    # Construct and print the final formula string
    # V = -(μ_2 + μ_3) / (2e)
    formula = (f"{voltage_symbol} = -({mu_symbol}_{stage_index_2} + "
               f"{mu_symbol}_{stage_index_3}) / "
               f"({denominator_coefficient}*e)")

    print("The simple formula that best approximates the second plateau is:")
    print(formula)

# Execute the function to display the result
print_voltage_formula()