import sympy

# (a) Diversity Order Calculation

# Define complex symbols for the analysis
# e1 represents the error symbol (x1 - x1')
# e2 represents the error symbol (x2 - x2')
e1, e2 = sympy.symbols('e1 e2', complex=True)

print("Part (a): Determining the Diversity Order")
print("==================================================")
print("The diversity order of a space-time code is determined by the rank of the codeword difference matrix, Delta_S = S(x) - S(x').")
print("A code has full diversity (order N for N transmit antennas) if this matrix is full rank for any pair of distinct codewords.")
print("For a 2x2 square matrix, this is equivalent to its determinant being non-zero for any non-zero error vector e = (e1, e2).")
print("\n")


# --- Analysis for Code Sa ---
print("Analyzing Code Sa:")
# The difference matrix for Sa is formed by replacing x_i with e_i
delta_Sa = sympy.Matrix([[e1, e2], [e2, e1]])
# Calculate the determinant
det_Sa = delta_Sa.det()

print("Difference matrix Delta_Sa:")
sympy.pprint(delta_Sa)
print(f"Determinant of Delta_Sa: {det_Sa}")
print("This determinant equals zero if e1**2 = e2**2, which means e1 = e2 or e1 = -e2.")
print("Since we can choose distinct codewords such that e1 = e2 (e.g., e1=1, e2=1),")
print("the determinant can be zero for a non-zero error vector.")
print("Therefore, the minimum rank is 1.")
print("Diversity Order for Sa is 1.")
print("-" * 50)


# --- Analysis for Code Sb ---
print("Analyzing Code Sb:")
# The difference matrix for Sb is formed by replacing x_i with e_i and x_i* with e_i*
delta_Sb = sympy.Matrix([[e1, e2], [e2, sympy.conjugate(e1)]])
# Calculate the determinant
det_Sb = delta_Sb.det()

print("Difference matrix Delta_Sb:")
sympy.pprint(delta_Sb)
# The determinant is |e1|^2 - e2^2
# sympy represents |e1|^2 as e1*conjugate(e1)
print(f"Determinant of Delta_Sb: {det_Sb}")
print("This determinant equals zero if Abs(e1)**2 = e2**2.")
print("This condition can be met for a non-zero error vector. For example, if we choose error")
print("symbols e1=2 and e2=2 (which can be formed from QAM constellations),")
print("then Abs(2)**2 - 2**2 = 4 - 4 = 0.")
print("Therefore, the minimum rank is 1.")
print("Diversity Order for Sb is 1.")
print("-" * 50)


# --- Analysis for Code Sc ---
print("Analyzing Code Sc:")
# The difference matrix for Sc
delta_Sc = sympy.Matrix([[-sympy.conjugate(e1), e2], [-sympy.conjugate(e2), -e1]])
# Calculate the determinant
det_Sc = delta_Sc.det()

print("Difference matrix Delta_Sc:")
sympy.pprint(delta_Sc)
# The determinant is |e1|^2 + |e2|^2
# sympy represents this as e1*conjugate(e1) + e2*conjugate(e2)
print(f"Determinant of Delta_Sc: {det_Sc}")
print("This determinant is the sum of squared magnitudes of the error symbols, |e1|^2 + |e2|^2.")
print("This sum is zero if and only if e1 = 0 and e2 = 0, which corresponds to a zero error vector (i.e., non-distinct codewords).")
print("For any distinct codewords (non-zero error vector), the determinant is strictly positive.")
print("Therefore, the matrix is always full rank (rank 2).")
print("Diversity Order for Sc is 2 (Full Diversity).")
print("-" * 50)
print("\n")


# (b) Maximum Directivity

print("Part (b): Code with Maximum Directivity")
print("==================================================")
print("In the context of space-time codes, 'directivity' can be interpreted as the code's effectiveness in mitigating the effects of fading.")
print("This effectiveness is primarily characterized by the diversity order.")
print("A higher diversity order leads to a steeper slope of the bit-error-rate curve at high SNR, signifying better performance.")
print("\nComparison:")
print(" - Sa: Diversity Order = 1")
print(" - Sb: Diversity Order = 1")
print(" - Sc: Diversity Order = 2")
print("\nCode Sc has the highest diversity order. Therefore, it provides the maximum 'directivity' or robustness against channel fading.")
