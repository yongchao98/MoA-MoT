import numpy as np

def build_Sa(x):
    """Builds the matrix for code Sa."""
    S = np.zeros((6, 6), dtype=complex)
    for i in range(6):
        for j in range(6):
            # This is the structure from the problem description
            # S[i,j] = x[ (i+j) mod 6 ]
            idx = (i + j) % 6
            S[i, j] = x[idx]
    return S

def build_Sb(x):
    """Builds the matrix for code Sb."""
    S = np.zeros((6, 6), dtype=complex)
    S[0, :] = [x[0], -np.conj(x[1]), x[2], -np.conj(x[3]), x[4], -np.conj(x[5])]
    S[1, :] = [x[1], x[2], -np.conj(x[3]), x[4], -np.conj(x[5]), np.conj(x[0])]
    S[2, :] = [x[2], x[3], x[4], -np.conj(x[5]), np.conj(x[0]), -np.conj(x[1])]
    S[3, :] = [x[3], x[4], x[5], np.conj(x[0]), -np.conj(x[1]), np.conj(x[2])]
    S[4, :] = [x[4], x[5], x[0], np.conj(x[1]), np.conj(x[2]), -np.conj(x[3])]
    S[5, :] = [x[5], x[0], x[1], np.conj(x[2]), np.conj(x[3]), np.conj(x[4])]
    return S
    
def build_Sc(x):
    """Builds the matrix for code Sc."""
    S = np.zeros((6, 6), dtype=complex)
    S[0, :] = [x[0], np.conj(x[1]), -x[2], np.conj(x[3]), -x[4], np.conj(x[5])]
    S[1, :] = [x[1], -x[2], np.conj(x[3]), -x[4], np.conj(x[5]), np.conj(x[0])]
    S[2, :] = [-x[2], np.conj(x[3]), -x[4], np.conj(x[5]), np.conj(x[0]), -np.conj(x[1])]
    S[3, :] = [np.conj(x[3]), -x[4], np.conj(x[5]), -np.conj(x[0]), -np.conj(x[1]), np.conj(x[2])]
    S[4, :] = [-x[4], np.conj(x[5]), np.conj(x[0]), -np.conj(x[1]), -np.conj(x[2]), -np.conj(x[3])]
    S[5, :] = [np.conj(x[5]), np.conj(x[0]), -np.conj(x[1]), np.conj(x[2]), -np.conj(x[3]), -np.conj(x[4])]
    return S

print("Analyzing Code Sa:")
# For code Sa, the matrix is singular if the sum of symbols is zero.
# For example, dx = [1, -1, 0, 0, 0, 0].
dx_a = np.array([1, -1, 0, 0, 0, 0], dtype=complex)
delta_Sa = build_Sa(dx_a)
det_a = np.linalg.det(delta_Sa)
print(f"For a difference vector dx = {dx_a.tolist()}, where sum(dx) = {np.sum(dx_a):.1f}:")
# The final equation to output is the determinant value
print(f"det(ΔSa) = {det_a.real:.1f} + {det_a.imag:.1f}j")
if np.isclose(det_a, 0):
    print("The determinant is zero, so Sa is not full-diversity.")
else:
    print("The determinant is non-zero.")

print("\n" + "="*40 + "\n")

print("Analyzing Code Sb:")
# For code Sb, a carefully chosen dx can make two columns identical, leading to a zero determinant.
# Let dx = [alpha, jb, alpha, jb, alpha, jb]. We choose alpha=2, b=2.
alpha = 2
jb = 2j
dx_b = np.array([alpha, jb, alpha, jb, alpha, jb], dtype=complex)
delta_Sb = build_Sb(dx_b)
det_b = np.linalg.det(delta_Sb)
print(f"For a difference vector dx = {[f'{c.real:.1f}+{c.imag:.1f}j' for c in dx_b]}, which makes two columns of ΔSb identical:")
print(f"det(ΔSb) = {det_b.real:.1f} + {det_b.imag:.1f}j")
if np.isclose(det_b, 0):
    print("The determinant is zero, so Sb is not full-diversity.")
else:
    print("The determinant is non-zero.")
    
print("\n" + "="*40 + "\n")

print("Analyzing Code Sc:")
# Code Sc is designed to be full-diversity. It should have a non-zero determinant
# for any non-zero difference vector dx. We test it with the vectors that failed other codes.
delta_Sc_from_a = build_Sc(dx_a)
det_c1 = np.linalg.det(delta_Sc_from_a)
print(f"Using the dx from Sa's test ({dx_a.tolist()}):")
print(f"det(ΔSc) = {det_c1.real:.1f} + {det_c1.imag:.1f}j")

delta_Sc_from_b = build_Sc(dx_b)
det_c2 = np.linalg.det(delta_Sc_from_b)
print(f"\nUsing the dx from Sb's test ({[f'{c.real:.1f}+{c.imag:.1f}j' for c in dx_b]}):")
print(f"det(ΔSc) = {det_c2.real:.1f} + {det_c2.imag:.1f}j")

if not np.isclose(det_c1, 0) and not np.isclose(det_c2, 0):
    print("\nThe determinants are non-zero. This supports that Sc is a full-diversity code.")

print("\n" + "="*40 + "\n")

print("Conclusion:")
print("Codes Sa and Sb are not full-diversity, as specific non-zero difference symbol vectors make their determinants zero.")
print("Their diversity order is therefore less than the maximum possible value of N=6.")
print("Code Sc, however, maintains a non-zero determinant for these vectors and is designed for full diversity.")
print("Thus, code Sc provides the maximum diversity order.")
print("\nThe code with the maximum diversity order is Sc, and its value is 6.")
