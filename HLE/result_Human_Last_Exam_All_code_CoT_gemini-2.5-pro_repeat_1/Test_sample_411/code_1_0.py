def solve_stbc_analysis():
    """
    Analyzes the diversity order and directivity of three given space-time block codes.
    """
    print("This script analyzes the diversity order of three space-time block codes (STBCs).")
    print("The diversity order is a crucial metric for performance in fading channels and is determined by the rank criterion.")
    print("According to the rank criterion, the diversity order is the minimum rank of the codeword difference matrix, ΔS = S(x) - S(e), over all pairs of distinct transmitted symbol vectors x and e.")
    print("Let the symbol differences be Δx₁ and Δx₂. Since the symbol vectors are distinct, at least one of Δx₁ or Δx₂ is non-zero.")
    print("-" * 70)

    # Part (a): Diversity Order Calculation
    print("\n(a) What is the diversity order for each code S_a, S_b, and S_c?\n")

    # Analysis for S_a
    print("--- Analysis of Code S_a = [[x₁, x₂], [x₂, x₁]] ---")
    print("The difference matrix is ΔS_a = [[Δx₁, Δx₂], [Δx₂, Δx₁]].")
    print("The determinant of ΔS_a determines its rank. det(ΔS_a) = (Δx₁)² - (Δx₂)².")
    print("The final equation for the determinant is: det = (1)*(Δx₁)² + (-1)*(Δx₂)².")
    print("The numbers in this equation are 1 and -1.")
    print("The determinant is zero if Δx₁ = ±Δx₂. This condition is possible for symbols from an M-QAM constellation.")
    print("When the determinant is zero, the rank of ΔS_a is 1. Therefore, the minimum rank is 1.")
    print("Result: Diversity order for S_a = 1.\n")

    # Analysis for S_b
    print("--- Analysis of Code S_b = [[x₁, x₂], [x₂, x₁*]] ---")
    print("The difference matrix is ΔS_b = [[Δx₁, Δx₂], [Δx₂, (Δx₁)*]].")
    print("The determinant is det(ΔS_b) = Δx₁*(Δx₁)* - (Δx₂)² = |Δx₁|² - (Δx₂)². Note: * denotes complex conjugate.")
    print("The final equation for the determinant is: det = (1)*|Δx₁|² + (-1)*(Δx₂)².")
    print("The numbers in this equation are 1 and -1.")
    print("The determinant can be zero. For example, if Δx₁ and Δx₂ are real-valued (from M-PAM) and equal, then |Δx₁|² = (Δx₂)².")
    print("When the determinant is zero, the rank of ΔS_b is 1. Therefore, the minimum rank is 1.")
    print("Result: Diversity order for S_b = 1.\n")

    # Analysis for S_c
    print("--- Analysis of Code S_c = [[-x₁*, x₂], [-x₂*, -x₁]] ---")
    print("The difference matrix is ΔS_c = [[-(Δx₁)*, Δx₂], [-(Δx₂)*, -Δx₁]].")
    print("The determinant is det(ΔS_c) = (-(Δx₁)*)(-Δx₁) - (Δx₂)(-(Δx₂)*) = |Δx₁|² + |Δx₂|².")
    print("The final equation for the determinant is: det = (1)*|Δx₁|² + (1)*|Δx₂|².")
    print("The numbers in this equation are 1 and 1.")
    print("The determinant |Δx₁|² + |Δx₂|² is a sum of non-negative terms. It is zero only if both Δx₁=0 and Δx₂=0.")
    print("This is not possible for distinct codewords. Thus, the determinant is always positive.")
    print("Since the determinant is never zero, the matrix ΔS_c is always full rank (rank 2).")
    print("Result: Diversity order for S_c = 2.\n")
    
    print("Summary for part (a):")
    print("The diversity order for S_a is 1.")
    print("The diversity order for S_b is 1.")
    print("The diversity order for S_c is 2.")
    print("-" * 70)
    
    # Part (b): Maximum Directivity
    print("\n(b) Which code provides the maximum directivity?\n")
    print("In the context of space-time coding, 'directivity' is interpreted as the diversity gain, which describes the system's robustness to channel fading.")
    print("A higher diversity order provides a higher diversity gain.")
    print("The maximum possible diversity order for a system with N=2 transmit antennas and L=1 receive antenna is N*L = 2.")
    print("From our analysis in part (a), code S_c achieves this maximum diversity order of 2.")
    print("Therefore, S_c provides the maximum directivity (diversity gain).")
    print("-" * 70)

if __name__ == "__main__":
    solve_stbc_analysis()