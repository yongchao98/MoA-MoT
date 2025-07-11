import math

# 1. Given Parameters
Gk = 1000.0  # kN (Permanent vertical load)
Qk_v = 1500.0 # kN (Variable vertical load)
Qk_h = 300.0  # kN (Variable horizontal load)
h_load = 2.0  # m (Height of load application above ground)
Df = 0.75    # m (Depth of footing base)
h_base = 0.75 # m (Thickness of footing base)
gamma_concrete = 24.0  # kN/m^3
phi_k = 35.0  # degrees (Characteristic friction angle)
c_k = 0.0     # kPa (Characteristic cohesion)
gamma_soil = 20.0  # kN/m^3

# ULS Partial Factors for DA1-1
gamma_G = 1.35  # for permanent actions
gamma_Q = 1.50  # for variable actions

# Assumed footing dimension from iterative design process
B = 2.3  # m
L = B    # m

print(f"Plan: Calculate the required ULS design bearing resistance (q_Ed) for an assumed footing size B = {B} m.")
print("-" * 30)

# 2. Calculate Design Loads and Moments
# Self-weight of the footing base
SW_footing = B * L * h_base * gamma_concrete
print(f"Self-weight of footing (SW) = {B} m * {L} m * {h_base} m * {gamma_concrete} kN/m^3 = {SW_footing:.2f} kN")

# Design vertical load (Vd)
Vd = gamma_G * (Gk + SW_footing) + gamma_Q * Qk_v
print(f"Design Vertical Load (Vd) = {gamma_G} * ({Gk} kN + {SW_footing:.2f} kN) + {gamma_Q} * {Qk_v} kN = {Vd:.2f} kN")

# Design horizontal load (Hd)
Hd = gamma_Q * Qk_h
print(f"Design Horizontal Load (Hd) = {gamma_Q} * {Qk_h} kN = {Hd:.2f} kN")

# Design moment at the base of the footing (Md)
lever_arm = h_load + Df
Md = Hd * lever_arm
print(f"Design Moment at Base (Md) = {Hd:.2f} kN * ({h_load} m + {Df} m) = {Md:.2f} kNm")
print("-" * 30)

# 3. Calculate Eccentricity and Effective Dimensions
# Eccentricity (e)
e = Md / Vd
print(f"Eccentricity (e) = {Md:.2f} kNm / {Vd:.2f} kN = {e:.4f} m")

# Effective width (B')
B_prime = B - 2 * e
print(f"Effective Width (B') = {B} m - 2 * {e:.4f} m = {B_prime:.4f} m")
L_prime = L  # No eccentricity in L direction

# Effective Area (A')
A_prime = B_prime * L_prime
print(f"Effective Area (A') = {B_prime:.4f} m * {L_prime} m = {A_prime:.4f} m^2")
print("-" * 30)

# 4. Calculate Required ULS Design Bearing Resistance (q_Ed)
q_Ed = Vd / A_prime
print("The required ULS design bearing resistance is the design pressure exerted by the footing.")
print(f"q_Ed = Vd / A' = {Vd:.2f} kN / {A_prime:.4f} m^2")
print(f"q_Ed = {q_Ed:.2f} kN/m^2")

# Final Answer
print("\nFinal Answer:")
print(f"The required ULS design bearing resistance is {q_Ed:.1f} kN/m^2.")