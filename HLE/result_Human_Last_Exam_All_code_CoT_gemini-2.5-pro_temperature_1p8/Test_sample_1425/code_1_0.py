def display_partition_function():
    """
    This function prints the formula for the grand canonical partition function
    using the path integral formulation and explains its components.
    """

    # Explanation of the terms
    explanation = """
The grand canonical partition function (Z) is expressed as a functional integral over all possible field configurations. Each configuration is weighted by a factor of exp(-S_E), where S_E is the Euclidean action.

Here is an explanation of each symbol in the final equation:
Z: The grand canonical partition function.
∫ D[ψ*, ψ]: The functional (or path) integral over all configurations of the field ψ(x, τ) and its complex conjugate ψ*(x, τ).
S_E: The Euclidean action, which governs the dynamics of the system in imaginary time.
exp(...): The exponential function.
∫[0 to β] dτ: Integral over imaginary time 'τ' from 0 to β.
∫ d^d x: Integral over 'd'-dimensional space.
β: The inverse temperature, defined as 1 / (k_B * T), where k_B is the Boltzmann constant and T is the temperature.
ψ(x, τ): The field describing the particles (bosonic or fermionic). It depends on position 'x' and imaginary time 'τ'.
ψ*(x, τ): The complex conjugate of the field.
∂/∂τ: The partial derivative with respect to imaginary time.
H(ψ*, ψ): The Hamiltonian density of the system, containing kinetic and potential/interaction terms.
μ: The chemical potential, which controls the average number of particles in the system.
"""

    # Final formula for the partition function Z and action S_E
    final_equation = """
The final expression for the partition function Z is:

Z = ∫ D[ψ*, ψ] exp(-S_E[ψ*, ψ])

Where the Euclidean action S_E is given by:

S_E = ∫[0 to β] dτ ∫ d^d x { ψ*(x, τ) * (∂/∂τ)ψ(x, τ) + H(ψ*, ψ) - μψ*(x, τ)ψ(x, τ) }
"""

    print(explanation)
    print(final_equation)

# Execute the function to display the result
display_partition_function()