def print_grassmann_integral_normalization():
    """
    This function prints the fundamental normalization rule for the
    Berezin integral over a single Grassmann variable. This rule
    is key to maintaining consistency with the Pauli exclusion principle.
    """

    # Define the components of the normalization equation: ∫ dη η = 1
    integral_symbol = "\u222B"  # Unicode for the integral sign
    measure_symbol = "dη"
    variable_symbol = "η"
    equals_symbol = "="
    normalization_value = 1

    print("In the path integral formalism, the Pauli exclusion principle is encoded")
    print("by the property that Grassmann variables square to zero (η² = 0).")
    print("\nThe integral measure 'dη' is defined by a set of rules, not a single value.")
    print("The key normalization rule that defines the integral is:")
    
    # We print the full equation as requested, highlighting each part.
    # The instruction was: "output each number in the final equation!"
    # We will print the equation and then the numerical value itself.
    print(f"{integral_symbol} {measure_symbol} {variable_symbol} {equals_symbol} {normalization_value}")

    print("\nThe specific 'value' that results from this normalized integral is:")
    print(normalization_value)

if __name__ == "__main__":
    print_grassmann_integral_normalization()