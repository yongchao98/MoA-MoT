def simulate_berezin_integral(integrand):
    """
    This function simulates the result of a Berezin integral for a single
    Grassmann variable η, based on its defining rules.

    Args:
        integrand (str): A string representing the function to be integrated.
                         Supported values are "1", "η", and "η^2".

    Returns:
        int: The result of the symbolic integral.
    """
    if integrand == "1":
        # The integral of a constant over a Grassmann variable is 0.
        return 0
    elif integrand == "η":
        # The integral of the variable itself is 1 by definition (normalization).
        return 1
    elif integrand == "η^2":
        # This case demonstrates the Pauli exclusion principle.
        # Since η^2 = 0 for any Grassmann variable, the integrand is 0.
        return 0
    else:
        raise ValueError("This simulation only supports integrands '1', 'η', and 'η^2'.")

def main():
    """
    Prints the fundamental rules of Grassmann integration, which encode the
    Pauli exclusion principle.
    """
    print("The Berezin integration measure 'dη' for a Grassmann variable 'η' is defined by rules that uphold the Pauli exclusion principle (η^2 = 0).")
    print("-" * 80)

    # Rule 1: Integral of a constant (unoccupied state component)
    integrand_1 = "1"
    result_1 = simulate_berezin_integral(integrand_1)
    # The equation contains the numbers 1 and 0.
    print(f"Rule 1: The integral over a constant term is zero.")
    print(f"   ∫ dη ⋅ {integrand_1} = {result_1}\n")

    # Rule 2: Normalization (occupied state component)
    integrand_2 = "η"
    result_2 = simulate_berezin_integral(integrand_2)
    # The equation contains the number 1.
    print(f"Rule 2: The normalization for the integral is one.")
    print(f"   ∫ dη ⋅ {integrand_2} = {result_2}\n")

    # Consequence of the rules and nilpotency (Pauli principle)
    integrand_3 = "η^2"
    result_3 = simulate_berezin_integral(integrand_3)
    # The equation contains the numbers 2 and 0.
    print("Demonstration of the Pauli exclusion principle:")
    print(f"   The integral of a doubly-occupied state component (η^2) is zero.")
    print(f"   ∫ dη ⋅ {integrand_3} = {result_3}   (because {integrand_3} itself is {result_3})\n")
    print(f"The value from Rule 2 ({result_2}) is the essential normalization of the measure for an occupied state.")


if __name__ == "__main__":
    main()