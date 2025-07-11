import math

def get_critical_exponent_nu():
    """
    Provides the value of the critical exponent ν for the G₄ (Φ⁴) theory.

    This function interprets the "G₄-theoretical framework" as the standard Φ⁴
    (phi-four) model, which is in the Ising universality class. The value of ν
    is dependent on the spatial dimension 'd'. This script will provide the
    high-precision value for the most physically relevant case, d=3.
    """

    # A dictionary of established values for the critical exponent ν for the
    # Ising / Φ⁴ universality class in different spatial dimensions 'd'.
    # The value for d=3 is a high-precision result from modern theoretical
    # and numerical methods (e.g., conformal bootstrap).
    nu_values = {
        2: 1.0,           # Exact value for d=2
        3: 0.629971,      # High-precision value for d=3
        4: 0.5            # Mean-field value for d>=4 (upper critical dimension)
    }

    # As the dimension 'd' was not specified, we select the canonical and
    # most physically significant non-trivial case: d = 3.
    selected_dimension = 3
    nu = nu_values[selected_dimension]

    # The final equation can be represented as: ν(d) = value
    # We will print each component of this relationship.
    print(f"Interpreting 'G₄-theory' as the standard Φ⁴ (phi-four) model:")
    print(f"For spatial dimension d = {selected_dimension},")
    print(f"The precise value of the critical exponent ν is: {nu}")


if __name__ == '__main__':
    get_critical_exponent_nu()