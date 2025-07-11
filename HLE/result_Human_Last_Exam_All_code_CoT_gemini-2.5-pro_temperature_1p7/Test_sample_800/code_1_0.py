import math

# Step 0: Define all given parameters and constants.
# Electrical parameters
P_el = 1500.0  # W
U_el = 220.0   # V
# Geometric parameters
D_R = 0.05 # m (5 cm)
# Temperature parameters in Celsius
theta_prime = 20.0    # Inlet air temperature
theta_double_prime = 60.0 # Outlet air temperature
theta_D = 180.0       # Wire temperature

print("--- Step 1: Calculation of Electrical Resistivity ---")
# Part a: a = f(pi) for f(t) in sin(t) = integral(e^(t-tau)f(tau)d(tau))
# Solving via Laplace Transform gives f(t) = cos(t) - sin(t).
a = math.cos(math.pi) - math.sin(math.pi)
# Part b: b = lim_{n->inf} [n^2 * integral_0^1 (x^n(1-x)) dx]
# The integral is 1/((n+1)(n+2)). The limit of n^2/(n^2+3n+2) as n->inf is 1.
b = 1.0
# Part c: c = 1/48 * integral_0^1 (ln(x))^4 dx
# The integral evaluates to the Gamma function Gamma(5) = 4! = 24.
c = 24.0 / 48.0
# Combine a, b, and c to find s and the electrical resistivity rho_el.
s = a + b + c
rho_el = s * 1e-6 # Electrical resistivity in Ohm*m

print(f"a = f(pi) = cos(pi) - sin(pi) = {a:.1f}")
print(f"b = lim [n^2 * integral] = {b:.1f}")
print(f"c = (1/48) * Gamma(5) = {c:.1f}")
print(f"s = a + b + c = {a:.1f} + {b:.1f} + {c:.1f} = {s:.1f}")
print(f"Electrical resistivity rho_el = s * 10^-6 = {rho_el:.1e} Ohm*m\n")

print("--- Step 2: Calculation of Airflow Properties ---")
# Use properties at the average air temperature (40 C) as per the hint.
theta_air_avg = (theta_prime + theta_double_prime) / 2.0
lambda_air = 27.354e-3  # W/(m K)
nu_air = 17.23e-6       # m^2/s
rho_air = 1.1124        # kg/m^3
Pr_air = 0.7056         # Prandtl number
cp_air = 1007.1         # J/(kg K)
# P_el = m_dot * cp * delta_theta_air
delta_theta_air = theta_double_prime - theta_prime
m_dot = P_el / (cp_air * delta_theta_air)
# w = m_dot / (rho * A_tube)
A_tube = math.pi * (D_R / 2.0)**2
w = m_dot / (rho_air * A_tube)
print(f"Mass flow rate m_dot = {P_el} / ({cp_air} * {delta_theta_air}) = {m_dot:.4f} kg/s")
print(f"Air velocity w = {m_dot:.4f} / ({rho_air} * {A_tube:.5f}) = {w:.3f} m/s\n")

print("--- Step 3: Calculation of Wire Length L ---")
# Two main equations link the wire's length L and diameter d.
# Eq 1 (Electrical): d^2 / L = (4 * rho_el * P_el) / (pi * U_el^2). This is a constant C1.
C1 = (4.0 * rho_el * P_el) / (math.pi * U_el**2)
# Eq 2 (Heat Transfer): P_el = h * A_surf * delta_theta_ht, with h from the Nusselt correlation.
# Combining the two equations and solving for L gives the final formula.
delta_theta_ht = theta_D - theta_air_avg
# This combined formula can be expressed as L^(5/4) = P_el / (K * C1^0.25)
# where K contains the heat transfer parameters.
K = (0.664 * w**0.5 * nu_air**-0.5 * Pr_air**(1.0/3.0) * lambda_air * math.pi * delta_theta_ht)
L_pow_5_over_4 = P_el / (K * C1**0.25)
L = L_pow_5_over_4**(4.0/5.0)
L_rounded = round(L)

print("The system of equations is solved for L.")
print(f"Intermediate term L^(5/4) = {L_pow_5_over_4:.2f} m^(5/4)")
print(f"Resulting wire length L = {L:.2f} m\n")
print("--- Final Answer ---")
print(f"The required length of the heating wire L, rounded to the nearest integer, is {L_rounded} m.")
print("\nFinal equation with all numerical values:")
# Prints the final full equation for L, with all values substituted.
print(f"L = ({P_el} / (0.664 * ({w:.3f})^0.5 * ({nu_air:.2e})^-0.5 * ({Pr_air:.4f})^(1/3) * {lambda_air:.5f} * {math.pi:.4f} * {delta_theta_ht} * (({4} * {rho_el:.1e} * {P_el}) / ({math.pi:.4f} * {U_el**2:.0f}))^0.25))^(4/5) = {L_rounded}")
<<<40>>>