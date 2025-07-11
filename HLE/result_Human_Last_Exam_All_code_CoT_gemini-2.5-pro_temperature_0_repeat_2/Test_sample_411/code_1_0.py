import sympy

def analyze_stbc():
    """
    Analyzes the diversity order of three space-time block codes (STBCs)
    and determines which one provides the maximum diversity.
    """
    # Define symbolic variables for two distinct sets of complex symbols
    x1, x2, e1, e2 = sympy.symbols('x1 x2 e1 e2', complex=True)

    # Define the error symbols, d1 = x1 - e1, d2 = x2 - e2.
    # For distinct codewords, at least one of d1 or d2 is non-zero.
    d1 = sympy.Symbol('d1', complex=True, nonzero=True)
    d2 = sympy.Symbol('d2', complex=True, nonzero=True)
    d1_z = sympy.Symbol('d1', complex=True) # version that can be zero
    d2_z = sympy.Symbol('d2', complex=True) # version that can be zero


    print("--- (a) Diversity Order Analysis ---")

    # --- Analysis of Code Sa ---
    print("\n--- Analysis of Code Sa = [[x1, x2], [x2, x1]] ---")
    # Difference matrix Delta_Sa = [[d1, d2], [d2, d1]]
    Delta_Sa = sympy.Matrix([[d1_z, d2_z], [d2_z, d1_z]])
    det_Sa = Delta_Sa.det()
    
    print("The determinant of the difference matrix is:")
    final_eq_sa = sympy.Eq(sympy.Symbol("det(ΔSa)"), det_Sa)
    # Output each number/symbol in the final equation
    print(f"{final_eq_sa.lhs} = {final_eq_sa.rhs.args[1]}**{final_eq_sa.rhs.args[1].exp} - {final_eq_sa.rhs.args[0]}**{final_eq_sa.rhs.args[0].exp}")


    print("\nTo find the diversity order, we check if the determinant can be zero for non-identical codewords.")
    print("The determinant is zero if d1**2 = d2**2, which means d1 = d2 or d1 = -d2.")
    print("This condition can be satisfied for non-zero d1 and d2 (e.g., if x1-e1 = x2-e2 != 0).")
    print("Since the difference matrix can be rank-deficient (rank 1), the diversity order for Sa is 1.")

    # --- Analysis of Code Sb ---
    print("\n--- Analysis of Code Sb = [[x1, x2], [x2, x1*]] ---")
    # Difference matrix Delta_Sb = [[d1, d2], [d2, d1*]]
    Delta_Sb = sympy.Matrix([[d1_z, d2_z], [d2_z, sympy.conjugate(d1_z)]])
    det_Sb = Delta_Sb.det()

    print("The determinant of the difference matrix is:")
    # sympy represents d1*conjugate(d1) as Abs(d1)**2
    final_eq_sb = sympy.Eq(sympy.Symbol("det(ΔSb)"), det_Sb)
    print(f"{final_eq_sb.lhs} = Abs({final_eq_sb.rhs.args[1].args[0]})**2 - {final_eq_sb.rhs.args[0].args[0]}**2")

    print("\nThe determinant is zero if Abs(d1)**2 = d2**2.")
    print("This can be satisfied for non-zero d1 and d2 (e.g., if d2 is real and d2 = Abs(d1)).")
    print("For example, using 4-QAM symbols, if x1=1+j, e1=-1+j, x2=1+j, e2=-1+j, then d1=2 and d2=2.")
    print("In this case, Abs(d1)**2 = 4 and d2**2 = 4, making the determinant zero.")
    print("Since the difference matrix can be rank-deficient (rank 1), the diversity order for Sb is 1.")

    # --- Analysis of Code Sc ---
    print("\n--- Analysis of Code Sc = [[-x1*, x2], [-x2*, -x1]] ---")
    # Difference matrix Delta_Sc = [[-d1*, d2], [-d2*, -d1]]
    Delta_Sc = sympy.Matrix([[-sympy.conjugate(d1_z), d2_z], [-sympy.conjugate(d2_z), -d1_z]])
    det_Sc = Delta_Sc.det()

    print("The determinant of the difference matrix is:")
    final_eq_sc = sympy.Eq(sympy.Symbol("det(ΔSc)"), det_Sc)
    print(f"{final_eq_sc.lhs} = Abs({final_eq_sc.rhs.args[0].args[0]})**2 + Abs({final_eq_sc.rhs.args[1].args[0]})**2")

    print("\nThe determinant is the sum of two non-negative terms: Abs(d1)**2 + Abs(d2)**2.")
    print("This sum is zero if and only if both Abs(d1)**2 and Abs(d2)**2 are zero, which means d1=0 and d2=0.")
    print("This corresponds to the case of identical codewords (S=E), which is excluded from the analysis.")
    print("For any distinct codewords, the determinant is always non-zero, so the matrix is always full rank (rank 2).")
    print("Therefore, the diversity order for Sc is 2.")

    print("\n\n--- (b) Code with Maximum Directivity ---")
    print("Maximum directivity corresponds to the maximum diversity order.")
    print("The diversity orders are: Sa=1, Sb=1, Sc=2.")
    print("The maximum possible diversity order for a 2x1 system is 2.")
    print("Code Sc achieves this maximum diversity order.")
    print("\nConclusion: Code Sc provides the maximum directivity.")

if __name__ == '__main__':
    analyze_stbc()