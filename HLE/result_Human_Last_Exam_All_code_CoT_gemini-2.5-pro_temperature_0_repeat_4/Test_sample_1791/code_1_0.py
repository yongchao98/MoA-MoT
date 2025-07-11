def display_effective_interaction():
    """
    Prints the formula for the effective electron-electron interaction
    mediated by phonons.
    """

    # Define the components of the formula as strings
    coupling_term = "g**2"
    momentum_term = "|q|**2"
    mass_term = "m"
    frequency_term = "w_q**2"

    # Construct the full formula string
    # The instruction is to output each "number" (term) in the equation.
    # We will print the final equation and then its components.
    
    print("The effective electron-electron interaction V_eff for a given momentum transfer q is:")
    
    # Print the final equation by combining the component strings
    final_equation = f"V_eff(q) = - ({coupling_term} * {momentum_term}) / ({mass_term} * {frequency_term})"
    print(final_equation)
    
    print("\nHere is a breakdown of each term in the equation:")
    print(f"Left-hand side (Interaction Potential): V_eff(q)")
    print(f"Sign: - (The interaction is attractive)")
    print(f"Numerator Part 1 (Coupling Constant Squared): {coupling_term}")
    print(f"Numerator Part 2 (Momentum Transfer Squared): {momentum_term}")
    print(f"Denominator Part 1 (Electron Mass): {mass_term}")
    print(f"Denominator Part 2 (Phonon Frequency Squared): {frequency_term}")

if __name__ == '__main__':
    display_effective_interaction()