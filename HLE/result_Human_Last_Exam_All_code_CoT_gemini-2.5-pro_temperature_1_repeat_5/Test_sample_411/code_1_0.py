import cmath

# This script analyzes the diversity order of three space-time block codes.

print("--- Diversity Order Analysis ---")
print("\nThe diversity order of an STBC in a 2x1 MISO system is the minimum rank of the codeword difference matrix, ΔS.")
print("A 2x2 matrix has rank 1 if its determinant is 0, and rank 2 if its determinant is non-zero.")
print("Let Δx1 = x1 - x'1 and Δx2 = x2 - x'2. For distinct codewords, at least one of Δx1 or Δx2 is non-zero.\n")

# --- Part (a): Diversity Order for each code ---

# Analysis for Code Sa
print("--------------------------------")
print("Analysis for Code Sa")
print("--------------------------------")
print("For S_a = [[x1, x2], [x2, x1]], the difference matrix is ΔS_a = [[Δx1, Δx2], [Δx2, Δx1]].")
print("The determinant is det(ΔS_a) = Δx1 * Δx1 - Δx2 * Δx2 = Δx1^2 - Δx2^2.")
print("This determinant can be zero if Δx1^2 = Δx2^2. This is possible if Δx1 = Δx2.")
dx1_a, dx2_a = 2, 2
det_a = dx1_a**2 - dx2_a**2
print(f"Example: Let Δx1 = {dx1_a} and Δx2 = {dx2_a}. The determinant is {dx1_a}^2 - {dx2_a}^2 = {det_a}.")
print("Since the determinant can be zero for non-zero symbol differences, the minimum rank is 1.")
print(">>> Diversity Order for S_a = 1\n")


# Analysis for Code Sb
print("--------------------------------")
print("Analysis for Code Sb")
print("--------------------------------")
print("For S_b = [[x1, x2], [x2, x1*]], the difference matrix is ΔS_b = [[Δx1, Δx2], [Δx2, Δx1*]].")
print("The determinant is det(ΔS_b) = Δx1 * Δx1* - Δx2 * Δx2 = |Δx1|^2 - Δx2^2.")
print("This determinant can be zero if |Δx1|^2 = Δx2^2. This is possible if, for example, both Δx1 and Δx2 are real and equal.")
dx1_b, dx2_b = 2, 2
det_b = abs(dx1_b)**2 - dx2_b**2
print(f"Example: Let Δx1 = {dx1_b} and Δx2 = {dx2_b}. The determinant is |{dx1_b}|^2 - {dx2_b}^2 = {det_b}.")
print("Since the determinant can be zero for non-zero symbol differences, the minimum rank is 1.")
print(">>> Diversity Order for S_b = 1\n")


# Analysis for Code Sc
print("--------------------------------")
print("Analysis for Code Sc")
print("--------------------------------")
print("For S_c = [[-x1*, x2], [-x2*, -x1]], the difference matrix is ΔS_c = [[-Δx1*, Δx2], [-Δx2*, -Δx1]].")
print("The determinant is det(ΔS_c) = (-Δx1*) * (-Δx1) - (Δx2) * (-Δx2*) = |Δx1|^2 + |Δx2|^2.")
print("This determinant is the sum of two non-negative terms. It is zero only if both Δx1 and Δx2 are zero, which corresponds to a zero difference matrix.")
print("For any distinct codewords, the determinant is always positive.")
dx1_c, dx2_c = 1+1j, 1-1j
det_c = abs(dx1_c)**2 + abs(dx2_c)**2
print(f"Example: Let Δx1 = {dx1_c} and Δx2 = {dx2_c}. The determinant is |{dx1_c}|^2 + |{dx2_c}|^2 = {abs(dx1_c)**2:.1f} + {abs(dx2_c)**2:.1f} = {det_c:.1f}.")
print("Since the determinant is never zero for non-zero differences, the minimum rank is 2.")
print(">>> Diversity Order for S_c = 2\n")


# --- Part (b): Code with Maximum Directivity ---
print("================================")
print("Conclusion")
print("================================")
print("(a) The diversity orders are:")
print("    - S_a: 1")
print("    - S_b: 1")
print("    - S_c: 2\n")
print("(b) The code providing maximum 'directivity' (diversity order) is the one that achieves the")
print("maximum possible diversity of 2 for this system. That code is S_c.")
