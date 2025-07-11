import sympy as sp

# Define symbols
z, h, k, mu, drho, g, S0, es = sp.symbols('z h k mu Delta_rho g S_0 e_s')
C1, C2, C3, C4 = sp.symbols('C1 C2 C3 C4')
H = k * h

# General solution for f(z)
# f(z) is the z-dependent part of the stream function psi = f(z)*sin(kx)
c = sp.cosh(k * z)
s = sp.sinh(k * z)
f = C1 * c + C2 * s + C3 * z * c + C4 * z * s

# Derivatives of f(z)
df = sp.diff(f, z)
d2f = sp.diff(df, z)
d3f = sp.diff(d2f, z)

# Boundary Conditions (BCs) expressed in terms of f(z)
# u_x = df/dz * sin(kx), u_z = -k*f*cos(kx)
# sigma_xz = mu * (d2f + k**2 * f) * sin(kx)
# sigma_zz = (mu/k * d3f - 3*mu*k*df) * cos(kx)

# BC 1: u_x = 0 at z = 0  (zero horizontal velocity at the top)
# This means df/dz at z=0 must be 0
bc1 = df.subs(z, 0)
# print(f"BC1 equation: {bc1} = 0")

# BC 2: sigma_zz = drho * g * es * cos(kx) at z = 0 (Normal stress at the top)
# We use Delta_rho as requested in the final expression for chi
sigma_zz_expr = (mu/k * d3f - 3*mu*k*df)
bc2_lhs = sigma_zz_expr.subs(z, 0)
bc2_rhs = drho * g * es
# print(f"BC2 equation: {bc2_lhs} = {bc2_rhs}")

# BC 3: sigma_xz = S0 * sin(kx) at z = h (Basal shear stress at the bottom)
sigma_xz_expr = mu * (d2f + k**2 * f)
bc3_lhs = sigma_xz_expr.subs(z, h)
bc3_rhs = S0
# print(f"BC3 equation: {bc3_lhs} = {bc3_rhs}")

# BC 4: sigma_zz = 0 at z = h (Zero normal stress at the bottom)
bc4 = sigma_zz_expr.subs(z, h)
# print(f"BC4 equation: {bc4} = 0")

# We have a system of 4 linear equations for C1, C2, C3, C4.
# Let's form the equations.
eq1 = sp.Eq(bc1, 0)
eq2 = sp.Eq(bc2_lhs, bc2_rhs)
eq3 = sp.Eq(bc3_lhs, bc3_rhs)
eq4 = sp.Eq(bc4, 0)

# Solve the system for the constants. We are particularly interested in the relation
# between es and S0, so we can treat them as knowns for now.
solution = sp.solve([eq1, eq2, eq3, eq4], (C1, C2, C3, C4))

# Now find the expression for chi
# chi = es * drho * g / S0
# From eq2, we have es = bc2_lhs / (drho*g)
# So, we need to express bc2_lhs in terms of S0.
# The solution dictionary expresses the Cs in terms of S0 and es. We can substitute
# the solution for the Cs back into the expression for es.
# Let's rearrange eq2 to express es in terms of the constants:
es_expr = bc2_lhs / (drho * g)

# Substitute the solved constants into this expression
es_solved = es_expr.subs(solution)

# Now we have an equation for es in terms of S0.
# We form chi = es * drho * g / S0
# We can find the ratio es/S0 from our solved equation.
chi_expr = sp.simplify((es_solved / S0) * (drho * g))

# The expression will be in terms of k and h. Let's substitute H=k*h for a cleaner look.
ch = sp.cosh(H)
sh = sp.sinh(H)
chi_final = sp.simplify(chi_expr.subs(sp.cosh(k*h), ch).subs(sp.sinh(k*h), sh))

# Print the final result for chi
print("The stable amplitude factor chi is given by:")
print(f"chi = {chi_final}")

# For clarity in the final output format requested by the prompt,
# let's format the sympy expression into a readable string.
numerator, denominator = sp.fraction(chi_final)
print("\nFormatted for final answer:")
print(f"chi = ({sp.pretty(numerator, use_unicode=False)}) / ({sp.pretty(denominator, use_unicode=False)})")

chi_pretty_num = sp.printing.pretty(numerator, use_unicode=False)
chi_pretty_den = sp.printing.pretty(denominator, use_unicode=False)

# The output format requests <<<answer content>>> at the end
final_answer_string = f"({chi_pretty_num}) / ({chi_pretty_den})"
print(f"<<<chi = {final_answer_string}>>>")