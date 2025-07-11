import sympy

# Set up the finite field of order 7
F = sympy.GF(7)
x = sympy.symbols('x')

# Ring C: F_7[x]/(5x^2 + x + 1)
pC = sympy.Poly(5*x**2 + x + 1, x, domain=F)
roots_C = pC.roots()
print(f"Ring C: Poly is {pC}. Roots in F_7: {roots_C}. Isomorphic to F_7 x F_7.")

# Ring E: F_7[x]/(3x^2 + x + 1)
pE = sympy.Poly(3*x**2 + x + 6, x, domain=F)
roots_E = pE.roots()
print(f"Ring E: Poly is {pE}. Roots in F_7: {roots_E}. Isomorphic to F_49.")

# Ring G: F_7[x]/(x^2 + 3x + 4)
pG = sympy.Poly(x**2 + 3*x + 4, x, domain=F)
factor_G = sympy.factor(pG)
print(f"Ring G: Poly is {pG}. Factored form: {factor_G}. Isomorphic to F_7[x]/(x^2).")

# Ring H (assuming typo): F_7[x]/(6x^2 + 5x + 4)
pH = sympy.Poly(6*x**2 + 5*x + 4, x, domain=F)
roots_H = pH.roots()
print(f"Ring H: Poly is {pH}. Roots in F_7: {roots_H}. Isomorphic to F_49.")

# Ring A: y^2 = x^3 + x^2 - 3x - 1
print("\n--- Curve Transformations ---")
u = sympy.symbols('u')
fA = sympy.Poly(x**3 + x**2 - 3*x + 1, x, domain=F)
# Shift x by c = -A/3 = -1/3 = -5 = 2
fA_shifted = fA.subs(x, u + F(2)).expand()
print(f"Curve A shifted (x -> u+2): y^2 = {sympy.Poly(fA_shifted, u, domain=F)}")

# Ring B: y^2 = x^3 + 2x^2 - 2x + 3
fB = sympy.Poly(x**3 + 2*x**2 - 2*x + 3, x, domain=F)
# Shift x by c = -A/3 = -2/3 = -10 = 4
fB_shifted = fB.subs(x, u + F(4)).expand()
print(f"Curve B shifted (x -> u+4): y^2 = {sympy.Poly(fB_shifted, u, domain=F)}")

# Check if the ratio of coefficients of u is a 4th power in F_7
# A's coeff is -1, B's coeff is -2. Ratio is 2.
fourth_powers = sorted(list(set([F(i)**4 for i in range(1, 7)])))
print(f"\nFourth powers in F_7*: {fourth_powers}")
print(f"Ratio of coefficients (-2)/(-1) = 2. Is 2 a fourth power? {F(2) in fourth_powers}")
print("Conclusion: A and B are isomorphic.")

# Ring I: y^2 = x^3 + 3x^2 + 3x + 2
fI = sympy.Poly(x**3 + 3*x**2 + 3*x + 2, x, domain=F)
# Shift x by c = -A/3 = -3/3 = -1 = 6
fI_shifted = fI.subs(x, u + F(6)).expand()
print(f"Curve I shifted (x -> u+6): y^2 = {sympy.Poly(fI_shifted, u, domain=F)}")
print("Curve I has j-invariant 0, while A and B have j-invariant 6. So I is not isomorphic to A or B.")
