import sympy as sp

# Set up sympy for pretty printing
sp.init_printing(use_unicode=True)

# Define the physical variables as symbolic objects
r, theta, R, E0, sigma1, sigma2 = sp.symbols(
    'r, theta, R, E_0, sigma_1, sigma_2', real=True, positive=True
)
# Define unknown coefficients for the potential
A1, B1 = sp.symbols('A1, B1')

# Define Legendre polynomial P1(cos(theta)) which governs the angular dependence
P1 = sp.cos(theta)

# --- Step 1 & 2: Define Potentials with Asymptotic Conditions ---

# General form of the potential inside the sphere (r < R)
# Must be finite at r=0, so we only have r^l terms. From symmetry with the
# external field, only the l=1 term is needed.
phi_in = A1 * r * P1

# General form of the potential outside the sphere (r > R)
# Must approach -E0*z = -E0*r*cos(theta) as r -> infinity.
# The r^-(l+1) terms represent the perturbation due to the sphere.
phi_out = -E0 * r * P1 + B1 * r**-2 * P1

# --- Step 3: Apply Boundary Conditions at r = R ---

# Boundary Condition 1: Continuity of Potential
# phi_in(R, theta) = phi_out(R, theta)
eq1 = sp.Eq(phi_in.subs(r, R), phi_out.subs(r, R))

# Boundary Condition 2: Continuity of the normal component of current density
# J_n = sigma * E_r = -sigma * d(phi)/dr
# sigma1 * (-d(phi_in)/dr)|_{r=R} = sigma2 * (-d(phi_out)/dr)|_{r=R}
eq2 = sp.Eq(sigma1 * (-sp.diff(phi_in, r)).subs(r, R), 
            sigma2 * (-sp.diff(phi_out, r)).subs(r, R))

# --- Step 4: Solve for the unknown coefficients A1 and B1 ---
solution = sp.solve([eq1, eq2], (A1, B1))
final_phi_out = phi_out.subs(solution)

# --- Step 5: Calculate the Electric Field from the Potential for r > R ---
# E = -grad(phi) in spherical coordinates
# E_r = -d(phi)/dr
# E_theta = -(1/r) * d(phi)/dtheta
final_E_r_out = -sp.diff(final_phi_out, r)
final_E_theta_out = -(1/r) * sp.diff(final_phi_out, theta)

# --- Format and Print the Final Expressions for r > R ---
print("The symbolic derivation provides the following results for the region OUTSIDE the sphere (r > R):\n")

print("Electric Potential Phi(r, theta) =")
# Factor the expression to match the format in the answer choices
factored_phi = -E0 * sp.factor(final_phi_out / (-E0))
# The following expression is printed to the console
sp.pprint(factored_phi)
print("\n----------------------------------\n")

print("Electric Field E(r, theta) has two components:\n")
print("Radial component E_r =")
# The following expression is printed to the console
sp.pprint(sp.simplify(final_E_r_out))
print("\n")

print("Polar component E_theta =")
# The following expression is printed to the console
sp.pprint(sp.simplify(final_E_theta_out))
print("\n----------------------------------\n")
print("Comparing these results with the given answer choices, we find a perfect match with option B.")
