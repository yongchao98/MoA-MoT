import sympy
from sympy import symbols, Matrix, conjugate, expand, Abs

def analyze_stbc():
    """
    Analyzes the diversity order of three space-time block codes Sa, Sb, and Sc.
    """
    #
    # Introduction to the method
    #
    print("--- Diversity Order Analysis of Space-Time Codes ---")
    print("The diversity order is determined by the rank of the codeword difference matrix, Delta_S.")
    print("A code has full diversity if det(Delta_S) is non-zero for any two distinct codewords.")
    print("Let dx1 = x1 - x1' and dx2 = x2 - x2'. For distinct codewords, at least one of dx1 or dx2 is non-zero.")
    print("-" * 60)

    #
    # Define symbolic variables
    #
    # We define them as real to make the output of |dx1|^2 cleaner,
    # but the logic holds for complex symbols dx1, dx2. The key part is |z|^2 >= 0.
    dx1_re, dx1_im, dx2_re, dx2_im = symbols('dx1_re dx1_im dx2_re dx2_im', real=True)
    dx1 = dx1_re + sympy.I * dx1_im
    dx2 = dx2_re + sympy.I * dx2_im
    
    # Let's use simpler complex symbols for print display
    dx1_s, dx2_s = symbols('Δx_1 Δx_2', complex=True, nonzero=True)

    #
    # --- Analysis of Code Sa ---
    #
    print("\n(1) Analysis of Code Sa = [[x1, x2], [x2, x1]]")
    Delta_S_a = Matrix([[dx1_s, dx2_s], [dx2_s, dx1_s]])
    det_a = Delta_S_a.det()

    print("The difference matrix is:")
    sympy.pprint(Delta_S_a, use_unicode=True)
    print("\nIts determinant is det(ΔS_a) = (Δx_1)^2 - (Δx_2)^2")
    
    print("\nThis determinant can be zero even for distinct codewords (Δx_1 != 0 or Δx_2 != 0).")
    print("For example, if we choose Δx_1 = 1 and Δx_2 = 1:")
    print("det(ΔS_a) = (1)^2 - (1)^2 = 1 - 1 = 0")
    print("Since the determinant can be zero for non-zero error vectors, the matrix is not always full rank.")
    print(">>> Diversity order of Sa is 1.")
    print("-" * 60)

    #
    # --- Analysis of Code Sb ---
    #
    print("\n(2) Analysis of Code Sb = [[x1, x2], [x2, x1*]]")
    Delta_S_b = Matrix([[dx1_s, dx2_s], [dx2_s, conjugate(dx1_s)]])
    # The determinant is |dx1|^2 - (dx2)^2
    
    print("The difference matrix is:")
    sympy.pprint(Delta_S_b, use_unicode=True)
    print("\nIts determinant is det(ΔS_b) = |Δx_1|^2 - (Δx_2)^2")

    print("\nThis determinant can also be zero for distinct codewords.")
    print("This occurs if |Δx_1|^2 = (Δx_2)^2. This requires Δx_2 to be a real number and |Δx_1| = |Δx_2|.")
    print("For example, if we choose Δx_1 = 2 (which is real) and Δx_2 = 2:")
    # Using specific values
    val_dx1 = 2
    val_dx2 = 2
    det_val = Abs(val_dx1)**2 - val_dx2**2
    print(f"det(ΔS_b) = |{val_dx1}|^2 - ({val_dx2})^2 = {Abs(val_dx1)**2} - {val_dx2**2} = {det_val}")
    print("Since the determinant can be zero for non-zero error vectors, the matrix is not always full rank.")
    print(">>> Diversity order of Sb is 1.")
    print("-" * 60)

    #
    # --- Analysis of Code Sc ---
    #
    print("\n(3) Analysis of Code Sc = [[-x1*, x2], [-x2*, -x1]]")
    Delta_S_c = Matrix([[-conjugate(dx1_s), dx2_s], [-conjugate(dx2_s), -dx1_s]])
    det_c = expand(Delta_S_c.det()) # Simplifies to |dx1|^2 + |dx2|^2
    
    print("The difference matrix is:")
    sympy.pprint(Delta_S_c, use_unicode=True)
    # Re-writing the sympy expression to be more explicit for the user.
    # The output of det_c is Δx_1*conjugate(Δx_1) + Δx_2*conjugate(Δx_2)
    print("\nIts determinant is det(ΔS_c) = (-Δx_1)*(-conj(Δx_1)) - (Δx_2)*(-conj(Δx_2)) = |Δx_1|^2 + |Δx_2|^2")

    print("\n|Δx_1|^2 and |Δx_2|^2 are always non-negative.")
    print("Therefore, their sum is zero if and only if |Δx_1|^2 = 0 AND |Δx_2|^2 = 0.")
    print("This implies Δx_1 = 0 and Δx_2 = 0, which means the codewords were not distinct.")
    print("For any two distinct codewords, the determinant is strictly positive and thus non-zero.")
    print(">>> Diversity order of Sc is 2.")
    print("-" * 60)

    #
    # --- Final Conclusion ---
    #
    print("\n--- Final Results ---")
    print("(a) The diversity order for each code is:")
    print("    - Sa: 1")
    print("    - Sb: 1")
    print("    - Sc: 2")
    print("\n(b) Which code provides the maximum directivity?")
    print("The maximum possible diversity order for this 2x1 MISO system is 2.")
    print("Code Sc achieves this maximum diversity and is therefore the best among the three.")
    print("The term 'directivity' likely refers to 'diversity' in this context.")
    
    # Format the final answer as requested
    final_answer = "<<<(a) Diversity Orders: Sa = 1, Sb = 1, Sc = 2. (b) Code Sc provides the maximum diversity.>>>"
    print(f"\n{final_answer}")

if __name__ == '__main__':
    analyze_stbc()