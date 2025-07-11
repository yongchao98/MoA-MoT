import sys

def display_fermionic_partition_function_formula():
    """
    This function prints the formula for the fermionic partition function Z
    in the imaginary time path integral formalism. It uses Unicode characters
    for mathematical symbols to improve readability.
    """

    print("The formula for the fermionic partition function Z in the imaginary time representation is given by a path integral:")
    print("Z = ∫ Dψ̄ Dψ exp(-S[ψ̄, ψ])")

    print("\nWhere S[ψ̄, ψ] is the Euclidean action. For non-interacting fermions, the action is:")
    
    # We will print the components of the action formula piece by piece
    # to show each part of the final equation clearly.
    sys.stdout.write("S = ")
    sys.stdout.write("∫₀^β dτ ")
    sys.stdout.write("∫ d³x ")
    sys.stdout.write("ψ̄(x, τ) [")
    sys.stdout.write("ħ(∂/∂τ)")
    sys.stdout.write(" - ")
    sys.stdout.write("(ħ²/2m)∇²")
    sys.stdout.write(" - ")
    sys.stdout.write("μ")
    sys.stdout.write("] ψ(x, τ)")
    print("\n")

    print("The terms in the action are:")
    print("  - τ: Imaginary time, integrated from 0 to β.")
    print("  - β: Inverse temperature (1 / k_B T).")
    print("  - ψ, ψ̄: Anticommuting Grassmann fields representing the fermions.")
    print("  - ħ: Reduced Planck constant.")
    print("  - m: Mass of the fermion.")
    print("  - μ: Chemical potential.")
    print("  - ∇²: Laplacian operator.")

    print("\nThe path integral is taken over all fields that satisfy the anti-periodic boundary condition in imaginary time:")
    print("ψ(x, β) = -ψ(x, 0)")

display_fermionic_partition_function_formula()