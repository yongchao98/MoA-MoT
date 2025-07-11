def calculate_critical_exponent_nu():
    """
    Calculates the precise value of the critical exponent ν within a G₄-theoretical
    framework at or above its upper critical dimension (d >= 4).

    In this regime, the system is described by mean-field theory. The critical
    exponent ν, which describes the scaling of the correlation length, takes on
    a universal value.
    """
    # In mean-field theory, the critical exponent ν is given by the fraction 1/2.
    numerator = 1
    denominator = 2

    # Calculate the value of ν
    nu_value = numerator / denominator

    # Print the final equation showing all numbers, as requested.
    print(f"Within the G₄-theoretical framework for d>=4 (mean-field theory), the critical exponent ν is calculated as:")
    print(f"{numerator} / {denominator} = {nu_value}")

if __name__ == "__main__":
    calculate_critical_exponent_nu()