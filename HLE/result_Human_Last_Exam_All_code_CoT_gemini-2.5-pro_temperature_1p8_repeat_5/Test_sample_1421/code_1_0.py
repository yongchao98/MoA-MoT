def fermionic_partition_function_formula():
    """
    This function prints the formula for the fermionic partition function Z
    in the imaginary time representation.
    """

    formula = """
The formula for the fermionic partition function Z in the imaginary time representation is given by a functional integral over Grassmann fields ψ(x, τ) and ψ̄(x, τ):

Z = ∫ Dψ̄ Dψ * exp(-S_E[ψ̄, ψ])

Where:

- ∫ Dψ̄ Dψ: Represents the functional or path integral over all possible configurations of the anti-commuting Grassmann fields ψ̄ and ψ.

- S_E[ψ̄, ψ]: Is the Euclidean action for the fermionic system. For a system of non-relativistic fermions, it is given by:

  S_E = ∫[0 to β] dτ ∫ d³x  ψ̄(x, τ) * [∂/∂τ - (∇²/2m) - μ] * ψ(x, τ)

  Let's break down the terms in the action:
  - τ: Imaginary time, integrated from 0 to β.
  - β: Inverse temperature, β = 1 / (k_B * T).
  - x: Spatial coordinate, integrated over all space.
  - ψ(x, τ), ψ̄(x, τ): The fermionic fields, which are functions of space and imaginary time.
  - ∂/∂τ: The partial derivative with respect to imaginary time.
  - ∇²: The Laplacian operator, representing the kinetic energy.
  - m: The mass of the fermion.
  - μ: The chemical potential.
  
- Boundary Conditions: A crucial aspect for fermions is that the Grassmann fields must satisfy anti-periodic boundary conditions in imaginary time due to the trace operation in the definition of Z:
  
  ψ(x, 0) = -ψ(x, β)
  ψ̄(x, 0) = -ψ̄(x, β)
"""
    print(formula)

fermionic_partition_function_formula()

final_answer = "Z = ∫ Dψ̄ Dψ * exp( - ∫[0 to β] dτ ∫ d³x  ψ̄(x, τ) * [∂/∂τ - (∇²/2m) - μ] * ψ(x, τ) )"
print(f"\nFinal Compact Formula:\n{final_answer}")