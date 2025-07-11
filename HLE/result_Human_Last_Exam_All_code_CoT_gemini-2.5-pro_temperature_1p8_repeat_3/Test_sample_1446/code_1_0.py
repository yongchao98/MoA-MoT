def solve_critical_exponent():
    """
    Calculates and displays the mean-field value for the critical exponent ν.

    Within a Ginzburg-Landau φ⁴ theory, the critical exponents are given by
    mean-field theory for spatial dimensions d ≥ 4. In this regime, the
    correlation length exponent ν has a precise, universal value.
    """
    # Define the components of the fractional value for ν
    numerator = 1
    denominator = 2

    # Calculate the decimal value of the exponent
    nu_value = numerator / denominator

    # Output the explanation and the final equation, including all numbers
    print("In a G₄ (φ⁴) theoretical framework for dimensions d ≥ 4 (the mean-field regime):")
    print("The precise value of the critical exponent ν is given by the equation:")
    print(f"ν = {numerator} / {denominator} = {nu_value}")

solve_critical_exponent()