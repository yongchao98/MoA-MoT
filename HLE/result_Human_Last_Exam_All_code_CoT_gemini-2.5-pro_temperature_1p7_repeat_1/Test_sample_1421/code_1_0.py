def display_fermionic_partition_function_formula():
    """
    This function prints the formula for the fermionic partition function Z
    in the imaginary time path integral formalism, along with explanations of its components.
    """

    # Define the components of the formula using unicode for better readability
    equation_main = "Z = ∫ Dψ̄ Dψ  exp(-Sₑ[ψ̄, ψ])"
    action_definition = "Sₑ[ψ̄, ψ] = ∫₀^β dτ ∫ dᵈx  ψ̄(x,τ) [ħ∂/∂τ + H] ψ(x,τ)"
    boundary_condition = "Boundary Conditions: ψ(x, 0) = -ψ(x, β)"

    print("The formula for the fermionic partition function Z in the imaginary time path integral representation is:")
    print("-" * 80)
    print(f"  {equation_main}")
    print("-" * 80)

    print("\nWhere the main components are:\n")

    print(f"1.  Euclidean Action (Sₑ):")
    print(f"    The action for the system in imaginary time is given by:")
    print(f"    {action_definition}\n")

    print(f"2.  Functional Integral (∫ Dψ̄ Dψ):")
    print(f"    This is an integral over all possible configurations of the fermionic fields.")
    print(f"    For fermions, these fields (ψ, ψ̄) are anti-commuting Grassmann variables.\n")

    print(f"3.  Anti-Periodic Boundary Conditions:")
    print(f"    A crucial feature for fermions, the trace operation in the partition function")
    print(f"    imposes anti-periodic boundary conditions on the fields in the imaginary time direction:")
    print(f"    {boundary_condition}\n")

    print("Definitions of each symbol in the formula:")
    print("{:<4} : {}".format("Z", "The Grand Canonical Partition Function."))
    print("{:<4} : {}".format("ψ, ψ̄", "Independent fermionic Grassmann fields."))
    print("{:<4} : {}".format("Sₑ", "The Euclidean Action of the system."))
    print("{:<4} : {}".format("∫", "Represents integration (functional, temporal, and spatial)."))
    print("{:<4} : {}".format("exp", "The exponential function."))
    print("{:<4} : {}".format("ħ", "The reduced Planck's constant."))
    print("{:<4} : {}".format("∂/∂τ", "The partial derivative with respect to imaginary time."))
    print("{:<4} : {}".format("τ", "Imaginary time, which runs from 0 to β."))
    print("{:<4} : {}".format("β", "The inverse temperature (β = 1 / (k_B * T))."))
    print("{:<4} : {}".format("d", "The number of spatial dimensions."))
    print("{:<4} : {}".format("H", "The single-particle Hamiltonian operator (e.g., kinetic energy minus chemical potential)."))

# Execute the function to print the formula
display_fermionic_partition_function_formula()