import sys

def solve_critical_exponent():
    """
    Calculates the value of the critical exponent nu (ν) in the G₄-theory.

    The G₄-theoretical framework is a mean-field theory (specifically, the Landau-Ginzburg
    ϕ⁴ theory) for second-order phase transitions. This theory provides exact results
    for systems with spatial dimensions d greater than the upper critical dimension, d_c = 4.

    The critical exponent ν describes how the correlation length (ξ) diverges as the
    temperature approaches the critical temperature.

    In the context of this mean-field theory, the value for ν is a precise rational number.
    """

    # The value of the critical exponent ν in mean-field theory is 1/2.
    numerator = 1
    denominator = 2

    # Calculate the value
    nu_value = numerator / denominator

    # Output the explanation and the equation with each number.
    print(f"In the G₄-theoretical framework (mean-field theory), the critical exponent ν is defined by the equation:")
    print(f"ν = {numerator} / {denominator}")
    print(f"The precise value of the critical exponent ν is: {nu_value}")

solve_critical_exponent()