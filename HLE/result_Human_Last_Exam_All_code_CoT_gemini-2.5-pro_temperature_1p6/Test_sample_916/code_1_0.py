import sympy as sp

# Define the symbolic variables used in the problem
R, N, N0, I0, i0, omega, t, g, mu0, alphaT, T, T0, Bs = sp.symbols(
    'R N N_0 I_0 i_0 omega t g mu_0 alpha_T T T_0 B_s'
)

# Step 1: Define the instantaneous current in the second coil
i_t = i0 * sp.sin(omega * t)
print("The instantaneous current in the second coil is i(t):")
sp.pprint(i_t)
print("\n" + "="*50 + "\n")

# Step 2: Define the magnetic field intensity H0 from the DC coil
# This is derived from Ampere's law for the given geometry.
H0 = (N0 * I0) / g
print("The magnetic field intensity from the DC coil is H0:")
sp.pprint(H0)
print("\n" + "="*50 + "\n")

# Step 3: Model the effective permeability, including all non-linear effects.
# We first define the temperature-corrected permeability, mu_T.
temp_correction_factor = (1 - alphaT * (T - T0))
mu_T = mu0 * temp_correction_factor
print("The permeability corrected for temperature is mu_T:")
sp.pprint(mu_T)
print("\n" + "="*50 + "\n")

# For the saturation effect, we approximate the field B in the formula with the ideal
# unsaturated field, B_ideal = mu0 * H0.
B_ideal = mu0 * H0
saturation_denominator = 1 + B_ideal / Bs
print("The denominator for the saturation effect is (1 + B_ideal / B_s):")
sp.pprint(saturation_denominator)
print("\n" + "="*50 + "\n")

# The total effective permeability is the combination of these two effects.
mu_eff = mu_T / saturation_denominator
print("The combined effective permeability is mu_eff:")
sp.pprint(mu_eff)
print("\n" + "="*50 + "\n")

# Step 4: Calculate the actual magnetic flux density B0 in the gap.
B0 = mu_eff * H0
print("The resulting magnetic flux density from the DC coil is B0 = mu_eff * H0:")
sp.pprint(B0)
print("\n" + "="*50 + "\n")

# Step 5: Calculate the instantaneous Lorentz force f_x(t) on the second coil.
# Force = N * (current) * (length of one turn) * B_0
# The length of one turn of the outer coil is 2*pi*R.
# A negative sign is introduced by convention.
fx_t = -N * i_t * (2 * sp.pi * R) * B0

print("The final equation for the instantaneous force f_x(t) is constructed from the components above.")
print("The equation is: f_x(t) = - (N) * (i(t)) * (2 * pi * R) * (B0)")
print("Substituting all the derived terms and simplifying, we get the final expression:\n")

# We can simplify the expression to make it clearer and match the format of the choices.
# Numerator of the final expression
final_numerator = 2 * sp.pi * R * N * mu0 * temp_correction_factor * N0 * I0 * i0 * sp.sin(omega * t)
# Denominator of the final expression
final_denominator = g**2 * saturation_denominator

final_expression = -final_numerator / final_denominator

# Pretty print the final expression
sp.pprint(final_expression, use_unicode=True)