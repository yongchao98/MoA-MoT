def solve_stbc_diversity():
    """
    Analyzes the diversity order of three space-time block codes
    and determines which one provides the maximum diversity.
    """
    print("(a) Analysis of Diversity Order")
    print("="*40)
    
    # --- Analysis for Code S_a ---
    print("Code S_a:")
    print("The space-time code is S_a = [[x₁, x₂], [x₂, x₁]].")
    print("The difference matrix is ΔS_a = [[δ₁, δ₂], [δ₂, δ₁]], where δ₁ = x₁ - x'₁ and δ₂ = x₂ - x'₂.")
    
    # Printing the final equation for the determinant
    print("The determinant of the difference matrix is:")
    print("det(ΔS_a) = (δ₁)*(δ₁) - (δ₂)*(δ₂)")
    print("          = δ₁² - δ₂²")
    
    print("\nAnalysis:")
    print("The determinant det(ΔS_a) is zero if δ₁² = δ₂², which means δ₁ = ±δ₂.")
    print("For M-QAM constellations, it is possible to find distinct symbols such that their differences satisfy this condition (e.g., δ₁ = δ₂ ≠ 0).")
    print("Since the determinant can be zero for non-zero error vectors, the matrix is not always full rank.")
    print("Conclusion: The diversity order for S_a is 1.")
    print("-" * 40)

    # --- Analysis for Code S_b ---
    print("Code S_b:")
    print("The space-time code is S_b = [[x₁, x₂], [x₂, x₁*]].")
    print("The difference matrix is ΔS_b = [[δ₁, δ₂], [δ₂, δ₁*]], where δ₁* is the complex conjugate of δ₁.")

    # Printing the final equation for the determinant
    print("The determinant of the difference matrix is:")
    print("det(ΔS_b) = (δ₁)*(δ₁*) - (δ₂)*(δ₂)")
    print("          = |δ₁|² - δ₂²")

    print("\nAnalysis:")
    print("The determinant det(ΔS_b) is zero if |δ₁|² = δ₂².")
    print("This condition can be met for non-zero error vectors. For example, if symbols are chosen from a BPSK constellation,")
    print("δ₁ and δ₂ can be real numbers. If we choose δ₁ = δ₂ ≠ 0, the determinant is zero.")
    print("Therefore, the matrix is not always full rank.")
    print("Conclusion: The diversity order for S_b is 1.")
    print("-" * 40)

    # --- Analysis for Code S_c ---
    print("Code S_c:")
    print("The space-time code is S_c = [[-x₁*, x₂], [-x₂*, -x₁]].")
    print("The difference matrix is ΔS_c = [[-δ₁*, δ₂], [-δ₂*, -δ₁]].")
    
    # Printing the final equation for the determinant
    print("The determinant of the difference matrix is:")
    print("det(ΔS_c) = (-δ₁*)*(-δ₁) - (δ₂)*(-δ₂*)")
    print("          = |δ₁|² + |δ₂|²")

    print("\nAnalysis:")
    print("The determinant det(ΔS_c) is the sum of two non-negative terms, |δ₁|² and |δ₂|².")
    print("This sum is zero if and only if |δ₁|² = 0 and |δ₂|² = 0, which means δ₁ = 0 and δ₂ = 0.")
    print("This corresponds to the zero error case (S₁ = S₂). For any two distinct codewords, the determinant is strictly positive.")
    print("Therefore, the difference matrix is always full rank.")
    print("Conclusion: The diversity order for S_c is 2.")
    print("=" * 40)
    
    print("\n(b) Code with Maximum Diversity")
    print("="*40)
    print("Comparing the diversity orders:")
    print("Diversity(S_a) = 1")
    print("Diversity(S_b) = 1")
    print("Diversity(S_c) = 2")
    print("\nThe maximum possible diversity for a 2x2 system is 2.")
    print("Code S_c is the only one that achieves the maximum diversity order.")

if __name__ == '__main__':
    solve_stbc_diversity()