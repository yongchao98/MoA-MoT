def print_force_law():
    """
    This function prints the derived force law for a thermally isolated
    freely jointed chain polymer.
    """

    print("The force law between the ends of the thermally isolated polymer is:")

    # Define the symbols as strings for readable output
    force_symbol = "F(x)"
    energy_symbol = "E(0)"
    extension_symbol = "x"
    segments_symbol = "n"
    length_symbol = "ell"

    # Define the numerical coefficients and powers from the derivation
    numerator_coeff = -2
    denominator_coeff_n = 1  # For n, the coefficient is 1
    denominator_coeff_ell = 1 # For ell, the coefficient is 1
    power_of_n = 2
    power_of_ell = 2

    # Print the equation, showing each numerical component explicitly as requested.
    # The negative sign indicates a restoring (attractive) force.
    print(f"{force_symbol} = ({numerator_coeff} * {energy_symbol} * {extension_symbol}) / ({segments_symbol}**{power_of_n} * {length_symbol}**{power_of_ell})")


if __name__ == "__main__":
    print_force_law()
