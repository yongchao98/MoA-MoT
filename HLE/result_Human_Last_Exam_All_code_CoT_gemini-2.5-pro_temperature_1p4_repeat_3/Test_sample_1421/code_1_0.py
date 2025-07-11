import sys

def display_fermionic_partition_function():
    """
    This function prints the formula for the fermionic partition function (Z)
    in the imaginary time path integral formalism.
    """

    # Using unicode characters for mathematical symbols for better representation.
    # Z = Partition Function
    # ∫ = Integral sign
    # Dψ Dψ̄ = Functional integration measure over Grassmann fields
    # exp = Exponential function
    # β = 1 / (k_B * T), inverse temperature
    # τ = Imaginary time variable
    # x = Spatial variable
    # ψ(x, τ), ψ̄(x, τ) = Fermionic Grassmann fields
    # ∂_τ = Partial derivative with respect to imaginary time
    # H₀ = Single-particle Hamiltonian (e.g., -ħ²∇²/2m)
    # μ = Chemical potential

    formula = (
        "\nThe formula for the fermionic partition function Z is given by the functional integral:\n\n"
        "         ╭             β           ╮\n"
        "Z = ∫ Dψ̄Dψ exp│- ∫ dτ ∫ d³x  L_E(ψ̄, ψ)│\n"
        "         ╰             0           ╯\n"
    )

    lagrangian = (
        "\nWhere the Euclidean Lagrangian density, L_E, is:\n\n"
        "L_E(ψ̄, ψ) = ψ̄(x, τ) (∂_τ + H₀ - μ) ψ(x, τ)\n"
    )

    full_formula = (
        "\nCombining these gives the full expression:\n\n"
        "          ╭β           ╮\n"
        "Z = ∫ Dψ̄Dψ exp │-∫dτ ∫d³x  ψ̄(x, τ)(∂_τ + H₀ - μ)ψ(x, τ) │\n"
        "          ╰0           ╯\n"
    )

    explanation = (
        "\nExplanation of the terms:\n"
        "  - Z: The partition function.\n"
        "  - ∫ Dψ̄Dψ: The functional integral over all configurations of the anticommuting Grassmann fields ψ̄ and ψ.\n"
        "  - β: Inverse temperature (β = 1/(k_B*T)). The time integral runs from 0 to β.\n"
        "  - τ: Imaginary time.\n"
        "  - ψ(x, τ), ψ̄(x, τ): The fermionic fields.\n"
        "  - ∂_τ: The partial derivative with respect to imaginary time.\n"
        "  - H₀: The single-particle Hamiltonian operator.\n"
        "  - μ: The chemical potential.\n"
    )

    boundary_conditions = (
        "\nCrucially, the integral is performed over fields that satisfy anti-periodic boundary conditions in imaginary time,\n"
        "which is a direct consequence of their fermionic nature:\n\n"
        "ψ(x, β) = -ψ(x, 0)\n"
    )

    print(formula)
    print(lagrangian)
    print(full_formula)
    print(explanation)
    print(boundary_conditions)

if __name__ == "__main__":
    display_fermionic_partition_function()
