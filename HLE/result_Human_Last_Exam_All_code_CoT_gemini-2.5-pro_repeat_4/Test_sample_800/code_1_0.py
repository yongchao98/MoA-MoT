import numpy as np

# Step 1: Define given constants and material properties
P_el = 1500.0  # W
U_el = 220.0  # V
D_R = 0.05  # m (Hairdryer tube diameter)
vartheta_in = 20.0  # °C
vartheta_out = 60.0  # °C
vartheta_D = 180.0  # °C (Wire temperature)

# Air properties at average temperature (40 °C)
lambda_m = 27.354e-3  # W/(m K)
nu_m = 17.23e-6  # m^2/s
rho_m = 1.1124  # kg/m^3
Pr_m = 0.7056
c_p_m = 1007.1  # J/(kg K)

# Step 2: Calculate the parameter 's'
# a = f(pi) where f(t) = cos(t) - sin(t)
a = np.cos(np.pi) - np.sin(np.pi)
# b = lim n->inf [n^2 * integral_0^1 x^n(1-x)dx] = lim n->inf [n^2 / ((n+1)(n+2))] = 1
b = 1.0
# c = 1/48 * integral_0^1 (ln(x))^4 dx = 1/48 * Gamma(5) = 1/48 * 4! = 24/48
c = 24.0 / 48.0
s = a + b + c
print(f"1. Calculating parameter s:")
print(f"   a = cos(π) - sin(π) = {a:.1f}")
print(f"   b = 1")
print(f"   c = 24 / 48 = {c:.1f}")
print(f"   s = a + b + c = {s:.1f}\n")

# Step 3: Calculate electrical properties of the wire
rho_el = s * 1e-6  # Ohm * m
R_el = U_el**2 / P_el  # Ohm
print(f"2. Calculating wire's electrical properties:")
print(f"   Electrical resistivity ρ_el = {s:.1f} * 10^-6 = {rho_el:.2e} Ω·m")
print(f"   Total resistance R_el = {U_el:.0f}² / {P_el:.0f} = {R_el:.3f} Ω\n")

# Step 4: Calculate air flow properties
vartheta_m = (vartheta_in + vartheta_out) / 2
delta_T_air = vartheta_out - vartheta_in
m_dot = P_el / (c_p_m * delta_T_air)
A_R = np.pi * D_R**2 / 4
w = m_dot / (rho_m * A_R)
print(f"3. Calculating air flow properties:")
print(f"   Average air temperature ϑ_m = ({vartheta_in:.0f} + {vartheta_out:.0f}) / 2 = {vartheta_m:.0f} °C")
print(f"   Air mass flow rate ṁ = {P_el:.0f} / ({c_p_m:.1f} * {delta_T_air:.0f}) = {m_dot:.4f} kg/s")
print(f"   Air velocity w = {m_dot:.4f} / ({rho_m:.4f} * {A_R:.5f}) = {w:.2f} m/s\n")

# Step 5: Formulate and solve the system of equations for L and d
# The two main equations relate L and d:
# Eq 1 (Resistance): L = (R_el * pi * d**2) / (4 * rho_el) => L = K2 * d**2
# Eq 2 (Heat Transfer): L = P_el / (h * pi * d * delta_T)
# where h = C_h * d**(-0.5). Substituting h gives: L = P_el / (C_h * pi * delta_T * d**0.5) => L = K1 / d**0.5

# Calculate the constants C_h, K1, and K2
delta_T_transfer = vartheta_D - vartheta_m
C_h = lambda_m * 0.664 * (w / nu_m)**0.5 * Pr_m**(1/3)
K1 = P_el / (C_h * np.pi * delta_T_transfer)
K2 = R_el * np.pi / (4 * rho_el)

# The solution to the system L = K2*d**2 and L = K1/d**0.5 is:
# L = K2**(1/5) * K1**(4/5)
L = K2**(1.0/5.0) * K1**(4.0/5.0)
L_rounded = round(L)

print(f"4. Solving for wire length L:")
print(f"   The system of equations for length L and diameter d is solved.")
print(f"   Calculated wire length L = ({K2:.2e})^(1/5) * ({K1:.4f})^(4/5)")
print(f"   L = {L:.3f} m\n")

print(f"5. Final Result:")
print(f"The required length of the heating wire is {L_rounded:.0f} m.")
print(f"Final Equation: L = (({R_el:.3f} * π) / (4 * {rho_el:.2e}))^(1/5) * ({P_el:.0f} / ({C_h:.3f} * π * {delta_T_transfer:.0f}))^(4/5) = {L_rounded:.0f}")

<<<15>>>