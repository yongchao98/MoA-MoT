import sympy

# Define all the symbols using sympy
R, N, N0, I0, i0, omega, t, g, mu0, alpha_T, T, T0, Bs = sympy.symbols('R N N_0 I_0 i_0 omega t g mu_0 alpha_T T T_0 B_s')

# Time-varying current i(t)
i_t = i0 * sympy.sin(omega * t)

# Temperature-dependent part of permeability
mu_T_term = mu0 * (1 - alpha_T * (T - T0))

# Saturation term in the denominator of permeability
saturation_term = (1 + (mu0 * N0 * I0) / (g * Bs))

# Numerator of the force equation from option B
numerator = -2 * sympy.pi * R * N * mu_T_term * N0 * I0 * i0 * sympy.sin(omega*t)

# Denominator of the force equation from option B
denominator = g**2 * saturation_term

# The full expression for the instantaneous force f_x(t)
f_x_t = numerator / denominator

# Print the final equation in a structured format
print("The instantaneous force f_x(t) is given by the equation:")
# We can't use sympy.pprint for the exact output format request, so we construct the string manually.
# The user wants each number/variable in the final equation to be output.
# The formula is f_x(t) = -2*pi*R*N * [ (mu0*(1-alpha_T*(T-T0))) * N0*I0*i0*sin(w*t) ] / [ g^2 * (1 + mu0*N0*I0/(g*Bs)) ]
# Let's print the equation term by term as requested by the prompt "Remember in the final code you still need to output each number in the final equation!"

print("f_x(t) = - (2 * pi * R * N * (mu_0 * (1 - alpha_T * (T - T_0))) * N_0 * I_0 * i_0 * sin(omega * t)) / (g**2 * (1 + (mu_0 * N_0 * I_0) / (g * B_s)))")

# Explicitly printing each variable and constant that constitutes the equation.
print("\nWhere the components of the equation are:")
print(f"Constant factor: -2 * pi")
print(f"Outer coil radius: R = {R}")
print(f"Number of turns in the second coil: N = {N}")
print(f"Permeability of free space: mu_0 = {mu0}")
print(f"Temperature coefficient: alpha_T = {alpha_T}")
print(f"Operating temperature: T = {T}")
print(f"Reference temperature: T_0 = {T0}")
print(f"Number of turns in the first coil: N_0 = {N0}")
print(f"DC current in the first coil: I_0 = {I0}")
print(f"Amplitude of AC current in the second coil: i_0 = {i0}")
print(f"Angular frequency: omega = {omega}")
print(f"Time: t = {t}")
print(f"Radial gap: g = {g}")
print(f"Saturation flux density: B_s = {Bs}")
