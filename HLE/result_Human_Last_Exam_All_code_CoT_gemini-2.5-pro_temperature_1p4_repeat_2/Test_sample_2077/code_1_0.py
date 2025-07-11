import sympy

# Define all the symbols
G_F, Q_W, F, M, E_nu, T, m_nu, pi = sympy.symbols('G_F Q_W F M E_nu T m_nu pi', positive=True, real=True)

# The differential cross section from Option F
# Denominator of the prefactor
den_F = pi * ((E_nu + M)**2 - (m_nu + M)**2) * ((E_nu + M)**2 - (m_nu - M)**2)
# Numerator of the prefactor
num_F = G_F**2 * Q_W**2 * F**2 * E_nu**2 * M**3
# The bracket term
bracket_F = 1 - T/E_nu - (M*T)/(2*E_nu**2) - m_nu**2/(2*E_nu**2) - (m_nu**2 * T)/(4*M*E_nu**2)

# Full expression for d(sigma)/dT from Option F
dsig_dT_F = (num_F / den_F) * bracket_F

print("Full expression from Option F:")
sympy.pprint(dsig_dT_F)
print("\n" + "="*50 + "\n")

# Step 1: Apply massless neutrino approximation (m_nu -> 0)
dsig_dT_m0 = dsig_dT_F.subs(m_nu, 0)
print("Expression after setting m_nu = 0:")
sympy.pprint(sympy.simplify(dsig_dT_m0))
print("\n" + "="*50 + "\n")


# Step 2: Apply low-energy approximation (E_nu << M)
# We can do this by taking the series expansion for E_nu -> 0
# Let's look at the prefactor and bracket separately.
prefactor_m0 = sympy.simplify(num_F / den_F.subs(m_nu, 0))
bracket_m0 = bracket_F.subs(m_nu, 0)

# The limit of the prefactor as E_nu -> 0
prefactor_approx = sympy.limit(prefactor_m0, E_nu, 0)

# For the bracket, we need to know the order of T. T is of order E_nu^2/M.
# So T/E_nu is of order E_nu/M, which goes to 0. MT/(2*E_nu^2) is of order 1.
# So we drop the T/E_nu term.
bracket_approx = 1 - (M*T)/(2*E_nu**2)

# The final approximated formula
dsig_dT_approx = prefactor_approx * bracket_approx

print("Final expression after E_nu << M approximation:")
sympy.pprint(dsig_dT_approx)
print("\n" + "="*50 + "\n")

# For comparison, let's write the original approximate formula from the problem
dsig_dT_original = (G_F**2 * M * Q_W**2 * F**2)/(4*pi) * (1 - (M*T)/(2*E_nu**2))
print("Original approximate formula from the problem description:")
sympy.pprint(dsig_dT_original)

# Verify they are the same
are_equal = sympy.simplify(dsig_dT_approx - dsig_dT_original) == 0
print(f"\nAre the derived approximation and the original formula equal? {are_equal}")
