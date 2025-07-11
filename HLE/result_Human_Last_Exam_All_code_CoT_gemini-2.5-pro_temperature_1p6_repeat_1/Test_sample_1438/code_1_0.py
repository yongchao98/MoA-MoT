def explain_grassmann_integral_measure():
    """
    Explains and prints the fundamental rules for Grassmann integration,
    which define the integral measure and enforce the Pauli exclusion principle.
    """

    # In fermionic path integrals, the Pauli exclusion principle is captured by
    # the property that Grassmann variables (eta) square to zero: eta^2 = 0.
    # The integration measure d_eta is defined by its results when integrating
    # the possible functions of eta. Since eta^2=0, any function of eta
    # can be written as F(eta) = a + b*eta.

    print("The integral measure for a Grassmann variable 'eta' is defined by the following rules (Berezin integration):")
    print("-" * 80)

    # Rule 1: The integral of a constant (or the unoccupied state) is defined to be zero.
    rule_1_value = 0
    print(f"Rule 1: The integral of 1 with respect to the measure d_eta is {rule_1_value}.")
    print(f"   Equation: ∫ d_eta = {rule_1_value}\n")

    # Rule 2: The integral of the variable itself (the occupied state) is defined to be one.
    # This provides the normalization for the measure.
    rule_2_value = 1
    print(f"Rule 2: The integral of eta with respect to the measure d_eta is {rule_2_value}.")
    print(f"   Equation: ∫ eta d_eta = {rule_2_value}\n")
    
    print("These rules mean the integral operation effectively extracts the coefficient of the 'eta' term.")
    print("For a function F(eta) = a + b*eta:")
    print("∫ (a + b*eta) d_eta = a * (∫ d_eta) + b * (∫ eta d_eta)")
    print(f"                      = a * ({rule_1_value}) + b * ({rule_2_value}) = b")

    # The question "what is the value of the measure" is best answered by the value
    # of the fundamental, non-trivial integral that normalizes it.
    print("-" * 80)
    print(f"The value of the fundamental integral ∫ eta d_eta, which ensures the Pauli exclusion principle, is {rule_2_value}.")

explain_grassmann_integral_measure()