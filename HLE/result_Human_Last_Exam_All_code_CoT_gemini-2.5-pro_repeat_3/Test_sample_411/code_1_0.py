def solve_stbc_diversity():
    """
    Analyzes and explains the diversity order for three given STBCs.
    The final conclusion is then printed.
    """
    print("This script analyzes the diversity order of three space-time block codes.")
    print("The diversity order is determined by the rank of the difference matrix ΔS = S(x1, x2) - S(x1', x2').")
    print("Full diversity of order 2 is achieved if det(ΔS) is non-zero for any non-zero error vector (Δx1, Δx2), where Δx1=x1-x1' and Δx2=x2-x2'.\n")

    # --- Analysis of Code Sa ---
    print("--- (a) Diversity Order Analysis ---\n")
    print("1. Code Sa = [[x1, x2], [x2, x1]]")
    print("   The difference matrix is ΔSa = [[Δx1, Δx2], [Δx2, Δx1]].")
    # Equation for the determinant
    print("   The determinant is calculated as: det(ΔSa) = (Δx1) * (Δx1) - (Δx2) * (Δx2)")
    print("   Final Equation: det(ΔSa) = (Δx1)^2 - (Δx2)^2")
    print("   This determinant is zero if Δx1 = Δx2 or Δx1 = -Δx2. This is possible for distinct codewords (e.g., if Δx1=1, Δx2=1).")
    diversity_a = 1
    print(f"   Therefore, the diversity order of Sa is {diversity_a}.\n")

    # --- Analysis of Code Sb ---
    print("2. Code Sb = [[x1, x2], [x2, x1*]]")
    print("   The difference matrix is ΔSb = [[Δx1, Δx2], [Δx2, (Δx1)*]].")
    # Equation for the determinant
    print("   The determinant is calculated as: det(ΔSb) = (Δx1) * (Δx1)* - (Δx2) * (Δx2)")
    print("   Final Equation: det(ΔSb) = |Δx1|^2 - (Δx2)^2")
    print("   This determinant can be zero if (Δx2)^2 is real and non-negative (i.e., Δx2 is real) and |Δx1| = |Δx2|.")
    print("   This condition can be met (e.g., for BPSK modulation, if Δx1=2, Δx2=2).")
    diversity_b = 1
    print(f"   Therefore, the diversity order of Sb is {diversity_b}.\n")

    # --- Analysis of Code Sc ---
    print("3. Code Sc = [[-x1*, x2], [-x2*, -x1]]")
    print("   The difference matrix is ΔSc = [[-(Δx1)*, Δx2], [-(Δx2)*, -Δx1]].")
    # Equation for the determinant
    print("   The determinant is calculated as: det(ΔSc) = (-(Δx1)*) * (-Δx1) - (Δx2) * (-(Δx2)*)")
    print("   Final Equation: det(ΔSc) = |Δx1|^2 + |Δx2|^2")
    print("   This is a sum of two non-negative terms. It can only be zero if both Δx1=0 and Δx2=0, which means the codewords are not distinct.")
    print("   Therefore, for any distinct pair of codewords, the determinant is non-zero.")
    diversity_c = 2
    print(f"   Therefore, the diversity order of Sc is {diversity_c} (full diversity).\n")

    # --- Answer for question (b) ---
    print("--- (b) Code with Maximum Diversity ---\n")
    print("The maximum possible diversity for a 2x1 MISO system is 2.")
    print(f"Comparing the diversity orders: Sa has {diversity_a}, Sb has {diversity_b}, and Sc has {diversity_c}.")
    print("Code Sc is the only one that achieves the maximum diversity order of 2.")

if __name__ == '__main__':
    solve_stbc_diversity()
