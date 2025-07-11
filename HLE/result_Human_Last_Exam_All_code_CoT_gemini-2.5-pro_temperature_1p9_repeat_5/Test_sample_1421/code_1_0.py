def print_fermionic_partition_function_formula():
    """
    Prints a textual representation of the formula for the fermionic partition function
    in the imaginary time path integral formalism.
    """
    # Define symbols for a cleaner presentation
    integral = "∫"
    psi = "ψ"
    psi_bar = "ψ̄"  # Represents the conjugate field
    beta = "β"      # Represents inverse temperature 1/(k_B*T)
    tau = "τ"       # Represents imaginary time
    partial = "∂"   # Represents the partial derivative symbol

    # Assemble the final formatted string containing the full formula
    formula_text = f"""
The formula for the fermionic partition function (Z) is given by the path integral over anticommuting Grassmann fields ({psi} and {psi_bar}):

Z = {integral} D[{psi_bar}]D[{psi}] exp(-S)

Where S is the Euclidean action. The action is an integral of the Euclidean Lagrangian (L_E) over imaginary time ({tau}) from 0 to {beta}, and over space (x):

               {beta}
              /
    S[{psi_bar}, {psi}] =   |  d{tau}  ∫ d³x  L_E
              /
              0

The Euclidean Lagrangian density, L_E, is given by:

    L_E = {psi_bar}(x, {tau}) [ħ * {partial}/{partial}{tau} + H] {psi}(x, {tau})

This path integral is evaluated subject to the crucial anti-periodic boundary conditions for fermions:

    {psi}(x, 0) = -{psi}(x, {beta})
    {psi_bar}(x, 0) = -{psi_bar}(x, {beta})

Where:
- H is the Hamiltonian of the system.
- ħ is the reduced Planck constant (often set to 1).
- {beta} is the inverse temperature (1 / k_B T).
"""

    print(formula_text)

# Execute the function to print the formula
print_fermionic_partition_function_formula()
