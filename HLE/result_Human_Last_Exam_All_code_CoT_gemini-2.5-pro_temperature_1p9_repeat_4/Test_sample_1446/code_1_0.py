import sys

def solve_critical_exponent():
    """
    This function explains and provides the precise value for the critical
    exponent nu (ν) in the context of a G₄-theoretical framework under
    mean-field theory conditions.
    """

    # The G₄-theoretical framework is typically a scalar φ⁴ field theory.
    # The critical exponents for this theory depend on the spatial dimension, d.
    # There exists an "upper critical dimension", d_c, at and above which
    # the exponents take on their simple mean-field values.
    # For a φ⁴ theory, the upper critical dimension is 4.
    d_upper_critical = 4

    # For any dimension d >= d_upper_critical, mean-field theory is exact.
    # We will state the precise value of ν from mean-field theory.
    # The exponent ν describes the divergence of the correlation length ξ:
    # ξ ~ |T - T_c|⁻ν

    # In mean-field theory, ν is precisely 1/2.
    numerator = 1
    denominator = 2
    nu_value = numerator / denominator

    print("For a G₄-theoretical framework (φ⁴ theory), the critical exponent ν depends on the spatial dimension d.")
    print(f"The upper critical dimension for this theory is d_c = {d_upper_critical}.")
    print("For d >= d_c, the system's behavior is described by mean-field theory, which provides a precise value for ν.")
    print("\nThe value is derived from the scaling of the correlation length, and the final equation for the exponent is:")

    # Outputting each number in the final equation as requested.
    print(f"ν = {numerator} / {denominator} = {nu_value}")

solve_critical_exponent()