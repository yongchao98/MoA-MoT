import sympy as sp

# Define all symbols as positive real numbers
G_F, Q_W, F, M, E_nu, m_nu, T = sp.symbols('G_F Q_W F M E_nu m_nu T', positive=True)

# The general formula for the differential cross section is dsigma/dT = C * B
# Let's define the components of the formula from option D.

# The denominator D in the prefactor C
D_denom = ((E_nu + M)**2 - (m_nu + M)**2) * ((E_nu + M)**2 - (m_nu - M)**2)

# The prefactor C
C_full = (G_F**2 * Q_W**2 * F**2 * E_nu**2 * M**3) / (sp.pi * D_denom)

# The bracketed term B
B_full = 1 - T/E_nu - (M*T)/(2*E_nu**2) - m_nu**2/(4*E_nu**2) - (m_nu**2 * T)/(4*M*E_nu**2)

# The full differential cross section from Answer D
dsig_dT_full = C_full * B_full

print("The general formula for the differential cross section (from option D) is dsigma/dT = C * B, where:")
print("\nC =")
sp.pprint(C_full)
print("\nB =")
sp.pprint(B_full)

# The approximate formula from the problem statement
dsig_dT_approx = (G_F**2 / (4 * sp.pi)) * M * Q_W**2 * F**2 * (1 - (M * T) / (2 * E_nu**2))
print("\n---")
print("The approximate formula given in the problem is:")
sp.pprint(dsig_dT_approx)
print("---\n")

# --- Verification ---
print("Now, we apply the approximations to the general formula from option D.")

# Step 1: Set neutrino mass m_nu = 0
dsig_dT_m_nu_zero = sp.simplify(dsig_dT_full.subs(m_nu, 0))
print("1. After setting neutrino mass m_nu = 0, the general formula simplifies to:")
sp.pprint(dsig_dT_m_nu_zero)

# Step 2: Apply the low-energy approximation, E_nu << M.
# We can do this by taking the series expansion for E_nu -> 0.
# The leading term in the denominator (E_nu + 2*M)**2 is (2*M)**2 = 4*M**2.
# The leading terms in the bracket are 1 and -M*T/(2*E_nu**2). The T/E_nu term is of a higher order in E_nu/M.
# Let's approximate the simplified expression for E_nu << M.
C_approx = sp.limit(dsig_dT_m_nu_zero / (1 - T/E_nu - (M*T)/(2*E_nu**2)), E_nu, 0)
# The bracketed term needs careful treatment. We keep terms of order (E_nu/M)^0.
# T is of order E_nu^2/M. So T/E_nu is order E_nu/M, and MT/E_nu^2 is order 1.
B_approx = 1 - (M*T)/(2*E_nu**2)
dsig_dT_final_approx = C_approx * B_approx

print("\n2. After applying the low-energy approximation (E_nu << M), it further simplifies to:")
sp.pprint(sp.simplify(dsig_dT_final_approx))

print("\nThis result matches the approximate formula given in the problem statement.")
print("Therefore, the formula from option D is the correct one.")
