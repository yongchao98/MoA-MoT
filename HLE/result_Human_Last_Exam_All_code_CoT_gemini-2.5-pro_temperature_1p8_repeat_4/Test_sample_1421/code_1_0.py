def print_fermionic_partition_function_formula():
    """
    Prints the formula for the fermionic partition function Z
    in the imaginary time path integral formalism.
    """
    
    # Define the components of the formula
    Z = "Z"
    integral_symbol = "∫"
    integration_measure = "Dψ̄Dψ"
    exponent = "e"
    action = "S_E[ψ̄, ψ]"
    action_definition_integral = "∫ dτ from 0 to β ∫ d³x"
    lagrangian_density = "L_E"
    lagrangian_expression = "ψ̄(τ,x)(∂_τ + H)ψ(τ,x)"
    boundary_condition = "ψ(τ, x) = -ψ(τ + β, x)"
    
    # Print the explanation and the components
    print("The formula for the fermionic partition function Z is a functional integral over Grassmann fields ψ and ψ̄.")
    print("The key components are:")
    print(f"Partition function: {Z}")
    print(f"Functional integral: {integral_symbol} {integration_measure}")
    print(f"Boltzmann factor: {exponent}^(-{action})")
    print("\nWhere the Euclidean action S_E is:")
    print(f"{action} = {action_definition_integral} {lagrangian_density}")
    print(f"And the Euclidean Lagrangian density (example for a single particle Hamiltonian H) is:")
    print(f"{lagrangian_density} = {lagrangian_expression}")
    print("\nThis integral is taken over fields that satisfy anti-periodic boundary conditions in imaginary time τ:")
    print(f"Boundary Condition: {boundary_condition}")
    
    # Assemble and print the final compact formula
    print("\n---")
    print("The final formula is:")
    # Print each part of the final equation
    print("Z = ", end="")
    print(f"{integral_symbol} ", end="")
    print(f"{integration_measure} ", end="")
    print("e^", end="")
    print("-[", end="")
    print("∫_0^β dτ ", end="")
    print("∫ d³x ", end="")
    print("ψ̄(∂_τ + H)ψ", end="")
    print("]")
    print("\n")


print_fermionic_partition_function_formula()

# Extracting the core formula for the final answer
# Z = ∫ Dψ̄Dψ exp[ -∫₀^β dτ ∫ d³x  ψ̄(τ,x)(∂_τ + H)ψ(τ,x) ]
# This must be accompanied by the anti-periodic boundary condition ψ(τ,x) = -ψ(τ+β, x)
final_formula = "Z = ∫ Dψ̄Dψ e^(-∫_0^β dτ ∫ d³x ψ̄(τ,x)(∂_τ + H)ψ(τ,x)) with boundary condition ψ(τ, x) = -ψ(τ + β, x)"
# Let's print the final answer tag below