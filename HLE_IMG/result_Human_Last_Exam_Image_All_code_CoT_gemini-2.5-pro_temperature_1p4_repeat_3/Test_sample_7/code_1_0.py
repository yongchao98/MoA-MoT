def print_voltage_formula():
    """
    This function derives and prints the simple formula for the second
    voltage plateau of a graphite anode based on the provided information.

    The second plateau corresponds to the phase transition between Stage 3 and Stage 2
    of lithium intercalation in graphite. The voltage of this plateau can be approximated
    by the average of the chemical potentials of the two coexisting stages.
    """

    # Define the numerical components of the formula.
    # The equation involves the chemical potentials of Stage 2 and Stage 3,
    # and division by 2 times the elementary charge 'e'. The leading
    # minus sign corresponds to a factor of -1.
    
    stage_index_1 = 2
    stage_index_2 = 3
    denominator_coefficient = 2
    
    # Define the symbolic components.
    voltage = "V"
    potential_symbol = "μ"
    charge_symbol = "e"
    
    # Construct and print the final formula.
    # The curly braces in the f-string will be replaced by the variable values.
    print("The formula for the second plateau is:")
    print(f"{voltage} = -({potential_symbol}_{stage_index_1} + {potential_symbol}_{stage_index_2}) / ({denominator_coefficient} * {charge_symbol})")

print_voltage_formula()