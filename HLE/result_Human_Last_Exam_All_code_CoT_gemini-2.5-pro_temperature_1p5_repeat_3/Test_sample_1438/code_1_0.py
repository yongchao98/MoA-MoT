def demonstrate_grassmann_integration_measure():
    """
    This function explains and demonstrates the defining properties of the
    Grassmann integral measure used in fermionic path integrals.
    """
    print("In fermionic path integrals, the Pauli exclusion principle is encoded by using anti-commuting Grassmann variables (η).")
    print("The key algebraic property of a Grassmann variable is:")
    
    # Representing the equation η * η = 0
    variable = "η"
    power = 2
    result_algebra = 0
    print(f"{variable}^{power} = {result_algebra}")
    
    print("\nThe integration over these variables (Berezin integration) is defined by a set of rules, which determines the 'measure' dη.")
    print("These rules are not based on a numerical value for dη, but on its operational definition.")
    print("The two fundamental rules for a single Grassmann variable are:\n")

    # Rule 1: The integral of a constant (or just the measure itself) is zero.
    integrand_1 = 1
    result_1 = 0
    print("Rule 1: The integral of the measure itself is 0.")
    print(f"∫ dη * {integrand_1} = {result_1}")
    
    print("-" * 25)

    # Rule 2: The integral is defined to be a projection onto the term linear in η.
    # This acts as a normalization condition.
    integrand_2 = "η"
    result_2 = 1
    print("Rule 2: The integral of the variable itself is 1.")
    print(f"∫ dη * {integrand_2} = {result_2}")
    
    print("\nThese two results, 0 and 1, are the values that define the integration measure.")

if __name__ == '__main__':
    demonstrate_grassmann_integration_measure()