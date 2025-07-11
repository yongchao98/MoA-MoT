def calculate_critical_exponent_nu():
    """
    Calculates and displays the critical exponent ν in the mean-field approximation.

    In the context of a G₄ (or φ⁴) theoretical framework, the value of the critical
    exponent ν depends on the spatial dimension 'd'. For dimensions d ≥ 4 (the upper
    critical dimension), mean-field theory provides the exact value. This script
    calculates and presents this classical value.
    """

    # In mean-field theory, the critical exponent ν is a simple rational number.
    # We define its numerator and denominator.
    numerator = 1
    denominator = 2

    # Calculate the precise value of ν.
    nu_value = numerator / denominator

    # Output the explanation and the final equation with its value.
    print("Within the mean-field approximation of the G₄-theoretical framework (valid for dimensions d >= 4):")
    print("The equation for the critical exponent ν is:")
    print(f"ν = {numerator} / {denominator}")
    print(f"The precise value of ν is: {nu_value}")

if __name__ == "__main__":
    calculate_critical_exponent_nu()