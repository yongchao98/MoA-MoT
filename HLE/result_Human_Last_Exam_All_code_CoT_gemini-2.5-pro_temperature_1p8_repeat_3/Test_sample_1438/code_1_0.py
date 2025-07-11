def explain_grassmann_integration():
    """
    Prints the defining rules of Berezin integration for a single Grassmann variable,
    which are essential for fermionic path integrals and the Pauli exclusion principle.
    """
    print("The Berezin integration measure for a Grassmann variable 'eta' is defined by two fundamental rules.")
    print("These rules are consistent with the Pauli exclusion principle (which implies eta^2 = 0).\n")

    # --- Rule 1: Integral of a constant ---
    # This corresponds to the contribution from an empty state.
    # The equation is: Integral(d_eta * 1) = 0
    integrand_1 = 1
    result_1 = 0
    print(f"Rule 1: The integral of a constant gives zero.")
    print(f"   Equation: Integral(d_eta * {integrand_1}) = {result_1}")
    print(f"   The numbers in this equation are: {integrand_1} and {result_1}\n")

    # --- Rule 2: Integral of the variable itself ---
    # This corresponds to the contribution from a singly-occupied state.
    # The equation is: Integral(d_eta * eta) = 1
    result_2 = 1
    print(f"Rule 2: The integral of the variable 'eta' itself is normalized to one.")
    print(f"   Equation: Integral(d_eta * eta) = {result_2}")
    print(f"   The number in this final equation is: {result_2}")
    print("\nThis value of 1 is the measure's normalization for an occupied state.")

explain_grassmann_integration()