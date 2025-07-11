def solve_electric_field():
    """
    This function prints the derived expressions for the electric field
    in the two regions of the cylindrical resistor.
    """
    # Symbolic representation of the variables
    V0 = "V_0"
    sigma1 = "sigma_1"
    sigma2 = "sigma_2"
    pi = "pi"
    r = "r"
    i_phi = "i_phi"

    # Construct the numerator and denominator for the E-field expressions
    # Numerator for E1 depends on sigma_2
    numerator_E1 = f"2 * {sigma2} * {V0}"

    # Numerator for E2 depends on sigma_1
    numerator_E2 = f"2 * {sigma1} * {V0}"

    # The denominator is the same for both regions
    denominator = f"r * {pi} * ({sigma1} + {sigma2})"

    # Print the final results
    print("Based on the derivation, the electric fields are:")
    print("-" * 50)
    
    print("In region 1 (0 < phi < pi/2):")
    # Print the equation for E1 with all its components
    print(f"E_1 = ({numerator_E1}) / ({denominator}) * {i_phi}")
    
    print("\nIn region 2 (pi/2 < phi < pi):")
    # Print the equation for E2 with all its components
    print(f"E_2 = ({numerator_E2}) / ({denominator}) * {i_phi}")
    print("-" * 50)

    print("\nThese expressions match answer choice C.")

if __name__ == "__main__":
    solve_electric_field()