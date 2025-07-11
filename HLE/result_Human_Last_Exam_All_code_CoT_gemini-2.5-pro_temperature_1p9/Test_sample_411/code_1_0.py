import sympy

# Set up symbolic variables
# e1 and e2 represent the difference between two distinct symbols (e.g., e1 = x1 - x1_hat).
# They are defined as complex since QAM symbols are complex numbers.
e1, e2 = sympy.symbols('e1 e2', complex=True)
e1_conj = sympy.conjugate(e1)
e2_conj = sympy.conjugate(e2)

print("Part (a): Determining the diversity order for each code.")
print("The diversity order is found by analyzing the determinant of the difference matrix C = S - S_hat.")
print("If det(C) is non-zero for all non-zero error vectors (e1, e2), the code has full diversity (order 2).\n")

# --- Analysis of Code Sa ---
# For Sa = [[x1, x2], [x2, x1]], the difference matrix Ca is [[e1, e2], [e2, e1]]
Ca = sympy.Matrix([[e1, e2], [e2, e1]])
det_Ca = sympy.det(Ca)
print("--- Analysis of Code Sa ---")
print("For Sa, the difference matrix is C_a = [[e1, e2], [e2, e1]].")
# The final equation for the determinant is: det(Ca) = 1*e1**2 - 1*e2**2
# Outputting each number in the final equation:
print("The determinant equation is: det(C_a) = (1)*(e1**2) + (-1)*(e2**2)")
print(f"This determinant, {det_Ca}, can be zero for non-zero errors (e.g., if e1 = 1 and e2 = 1).")
print("Since the determinant can be zero for distinct codewords, the code does not achieve full diversity.")
print("Diversity Order for Sa = 1\n")


# --- Analysis of Code Sb ---
# For Sb = [[x1, x2], [x2, x1*]], the difference matrix Cb is [[e1, e2], [e2, conjugate(e1)]]
Cb = sympy.Matrix([[e1, e2], [e2, e1_conj]])
det_Cb = sympy.det(Cb)
print("--- Analysis of Code Sb ---")
print("For Sb, the difference matrix is C_b = [[e1, e2], [e2, conjugate(e1)]].")
# The final equation for the determinant is: det(Cb) = 1*|e1|**2 - 1*e2**2
print("The determinant equation is: det(C_b) = (1)*(|e1|**2) + (-1)*(e2**2)")
print(f"This determinant, {det_Cb}, can also be zero for non-zero errors (e.g., if e1 = 1 and e2 = 1).")
print("Since the determinant can be zero for distinct codewords, this code also does not achieve full diversity.")
print("Diversity Order for Sb = 1\n")


# --- Analysis of Code Sc ---
# For Sc = [[-x1*, x2], [-x2*, -x1]], the difference matrix Cc is [[-conj(e1), e2], [-conj(e2), -e1]]
Cc = sympy.Matrix([[-e1_conj, e2], [-e2_conj, -e1]])
det_Cc = sympy.simplify(sympy.det(Cc)) # Simplifies to |e1|^2 + |e2|^2
print("--- Analysis of Code Sc ---")
print("For Sc, the difference matrix is C_c = [[-conjugate(e1), e2], [-conjugate(e2), -e1]].")
# The final equation for the determinant is: det(Cc) = 1*|e1|**2 + 1*|e2|**2
print("The determinant equation is: det(C_c) = (1)*(|e1|**2) + (1)*(|e2|**2)")
print(f"This determinant, {det_Cc}, is a sum of squared magnitudes.")
print("It is zero only if both |e1|^2=0 and |e2|^2=0, which implies e1=0 and e2=0.")
print("This means det(C_c) is non-zero for any pair of distinct codewords.")
print("Therefore, this code achieves full diversity.")
print("Diversity Order for Sc = 2\n")


# --- Conclusion ---
print("----------------------------")
print("           RESULTS          ")
print("----------------------------")
print("(a) Summary of Diversity Orders:")
print("    - Code Sa: 1")
print("    - Code Sb: 1")
print("    - Code Sc: 2")
print("\n(b) Code with Maximum Directivity:")
print("Assuming 'maximum directivity' refers to the highest diversity order, the best code is Sc.")
print("It achieves the maximum possible diversity of 2 for a 2x1 MISO system.")
