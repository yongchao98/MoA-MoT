def grassmann_integral_normalization():
    """
    This function demonstrates the fundamental normalization rule for Grassmann variable integration,
    which is essential for maintaining the Pauli exclusion principle in fermionic path integrals.

    The Pauli principle is encoded in the property η² = 0 for any Grassmann variable η.
    The integration over such variables is defined by the rules:
    1. ∫ dη = 0
    2. ∫ dη η = 1

    This second rule provides the normalization for the measure 'dη'. We will print this
    defining equation and its value.
    """

    # The value of the normalized integral of a Grassmann variable is 1 by definition.
    integral_value = 1
    
    # We represent the equation symbolically.
    # integral_symbol = "∫"
    # measure = "dη"
    # variable = "η"
    
    print("The fundamental normalization equation for the Grassmann integral measure is:")
    # The final equation is ∫ dη η = 1
    print(f"∫ dη η = {integral_value}")

grassmann_integral_normalization()