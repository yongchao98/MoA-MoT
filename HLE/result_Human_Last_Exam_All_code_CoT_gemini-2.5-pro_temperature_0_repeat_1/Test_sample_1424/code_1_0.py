def print_green_function_formula():
    """
    This function prints the formula for the bare Green's function,
    showing its dependence on single-particle energy eigenvalues.
    """

    # Define the symbols used in the equation as string variables
    # for clear presentation.
    green_function = "G_0(k, ω)"
    frequency = "ω"
    energy_eigenvalue = "ϵ_k"
    infinitesimal = "iδ"
    numerator = "1"

    # Explain the context and the components of the formula
    print("In the Feynman path integral formalism, the bare Green's function (G_0) describes the propagation of a non-interacting particle.")
    print("Its functional dependence on the single-particle energy eigenvalues (ϵ_k) in frequency-momentum space is given by the following relation:\n")

    # Print each component of the equation as requested
    print("--- Equation Components ---")
    print(f"Numerator: {numerator}")
    print(f"Frequency variable: {frequency}")
    print(f"Single-particle energy eigenvalue: {energy_eigenvalue}")
    print(f"Infinitesimal for causality: {infinitesimal}")
    print("---------------------------\n")

    # Construct and print the final equation
    denominator = f"({frequency} - {energy_eigenvalue} + {infinitesimal})"
    equation = f"{green_function} = {numerator} / {denominator}"

    print("The final equation is:")
    print(equation)

    # Provide an explanation of the terms
    print("\nWhere:")
    print(f" - {green_function}: The bare Green's function for a particle with momentum 'k' and frequency 'ω'.")
    print(f" - {frequency}: The energy (frequency) of the propagating particle.")
    print(f" - {energy_eigenvalue}: The energy of the single-particle eigenstate with momentum 'k'. This is the energy the particle would have in the absence of interactions.")
    print(f" - {infinitesimal}: An infinitesimal imaginary term that enforces the correct time-ordering (causality). Its sign depends on whether the energy level ϵ_k is occupied or empty in the ground state.")
    print("\nThis equation shows that the Green's function has a pole (diverges) when the particle's energy ω is equal to the eigenenergy ϵ_k, which corresponds to a stable particle excitation.")

# Execute the function to display the answer
print_green_function_formula()