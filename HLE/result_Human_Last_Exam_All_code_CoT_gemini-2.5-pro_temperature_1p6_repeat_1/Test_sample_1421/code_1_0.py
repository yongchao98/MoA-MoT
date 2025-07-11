def print_fermionic_partition_function_formula():
    """
    Prints the formula for the fermionic partition function Z
    in the imaginary time path integral formalism.
    """

    formula = """
The fermionic partition function Z in the imaginary time representation is given by a functional integral over anti-commuting Grassmann fields ψ*(τ) and ψ(τ):

  Z = ∫ D[ψ*, ψ] exp(-S[ψ*, ψ])

Where:
  - ∫ D[ψ*, ψ] represents the functional integral over all possible field configurations.
  - S[ψ*, ψ] is the Euclidean action.

The Euclidean action S is defined as:

  S[ψ*, ψ] = ∫[from 0 to β] dτ [ Σ_i ψ_i*(τ) (∂_τ)ψ_i(τ) + H(ψ*(τ), ψ(τ)) ]

Here:
  - τ is the imaginary time.
  - β = 1/(k_B T) is the inverse temperature. The integral runs from 0 to β.
  - ψ*(τ) and ψ(τ) are the Grassmann fields, representing the fermions.
  - ∂_τ is the partial derivative with respect to imaginary time.
  - H(ψ*, ψ) is the Hamiltonian of the system expressed in terms of the Grassmann fields.

A crucial element for fermions is the anti-periodic boundary condition imposed on the fields:

  Boundary Condition: ψ(τ) = -ψ(τ + β)
"""
    print(formula)

# Execute the function to print the formula
print_fermionic_partition_function_formula()