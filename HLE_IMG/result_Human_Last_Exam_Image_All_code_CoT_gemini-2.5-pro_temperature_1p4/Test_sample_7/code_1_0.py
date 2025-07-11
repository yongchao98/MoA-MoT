def print_voltage_formula():
    """
    Prints the formula for the second voltage plateau of graphite anode.
    
    The formula is derived based on the mean-field theory of staging transitions
    in graphite intercalation compounds.

    V_plateau2: The voltage of the second plateau (S3 -> S2 transition)
    mu_1: Chemical potential related to Stage 1
    mu_2: Chemical potential related to Stage 2
    e: The elementary charge
    """
    
    # The derived formula for the second plateau voltage (V_2) is:
    # V_2 = (mu_1 - 2*mu_2) / e
    # We will print this, explicitly showing the coefficients as requested.
    
    formula = "V_plateau2 = (1 * mu_1 - 2 * mu_2) / e"
    print(formula)

print_voltage_formula()