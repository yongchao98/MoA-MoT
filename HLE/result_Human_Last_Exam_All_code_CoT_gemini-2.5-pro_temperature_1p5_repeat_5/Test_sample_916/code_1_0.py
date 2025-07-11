import sympy as sp

# Define all the symbols
f_x, t, R, N, N0, I0, i0, omega, g, mu0, alpha_T, T, T0, Bs = sp.symbols(
    'f_x t R N N_0 I_0 i_0 omega g mu_0 alpha_T T_1 T_0 B_s'
)

# Current i(t)
i_t = i0 * sp.sin(omega * t)

# Temperature correction term
temp_term = (1 - alpha_T * (T - T0))

# Saturation term in the denominator
saturation_term = (1 + (mu0 * N0 * I0) / (g * Bs))

# Numerator of the main fraction
numerator = mu0 * temp_term * N0 * I0 * i_t

# Denominator of the main fraction
denominator = g**2 * saturation_term

# Complete expression for the force
force_expression = -2 * sp.pi * R * N * (numerator / denominator)

# Create the equation string
equation = sp.Eq(f_x, force_expression)

# Print the equation in a readable format
# We manually construct the string to match the desired output format exactly
# and print each component of the equation as requested.
print("The instantaneous force f_x(t) is given by the equation:")
print(f"f_x(t) = -2*\u03c0*R*N * ( (\u03bc_0 * (1 - \u03b1_T*(T_1 - T_0)) * N_0*I_0 * i_0*sin(\u03c9*t)) / (g^2 * (1 + (\u03bc_0*N_0*I_0)/(g*B_s))) )")
print("\nWhere each term represents:")
print("Force: f_x(t)")
print("Constant: -2*\u03c0")
print("Outer coil radius: R")
print("Inner coil turns: N")
print("Outer coil turns: N_0")
print("DC current: I_0")
print("AC current amplitude: i_0")
print("AC current waveform: sin(\u03c9*t)")
print("Permeability of free space: \u03bc_0")
print("Temperature coefficient: \u03b1_T")
print("Operating temperature: T_1")
print("Reference temperature: T_0")
print("Radial gap: g")
print("Saturation flux density: B_s")