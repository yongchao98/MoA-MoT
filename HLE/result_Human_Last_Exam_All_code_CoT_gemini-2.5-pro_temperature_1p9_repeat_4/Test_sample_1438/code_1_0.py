import sys

def display_berezin_integration_rules():
    """
    In the path integral formalism, the "measure" for integrating over a
    Grassmann variable is not a single number but is defined by a set of axiomatic rules.
    These rules, known as Berezin integration, are constructed to be consistent with
    the algebraic property η² = 0, which enforces the Pauli exclusion principle.

    This script prints the two fundamental rules that define this integral for a single
    Grassmann variable η.
    """

    print("The Berezin integral over a single Grassmann variable η is defined by the following two rules:")
    print("-" * 80)

    # Rule 1: The integral of a constant is zero.
    # This rule ensures that only terms containing the variable contribute.
    constant_part = 1
    result_1 = 0
    print(f"Rule 1: The integral of a constant ({constant_part}) is defined to be {result_1}.")
    print(f"∫ dη * {constant_part} = {result_1}")
    print("\n")


    # Rule 2: The integral of the variable itself is normalized to one.
    # This is the normalization convention.
    linear_part_coeff = 1
    result_2 = 1
    print(f"Rule 2: The integral of the variable η (multiplied by coefficient {linear_part_coeff}) is defined to be {result_2}.")
    print(f"∫ dη * η = {result_2}")
    print("-" * 80)
    print("\nThese two rules combined define the integral for any function f(η) = A + Bη as ∫ dη f(η) = B.")


if __name__ == "__main__":
    display_berezin_integration_rules()
