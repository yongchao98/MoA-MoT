import math

# Step 1: Define Inputs and Assumptions
# Given characteristic loads
Gk_structure = 1000  # kN (Permanent vertical load from superstructure)
Qk_v = 1500  # kN (Variable vertical load)
Qk_h = 300   # kN (Variable horizontal load)

# Given dimensions and material properties
h_column = 2.0  # m (Height of column above ground)
D_footing = 0.75  # m (Depth of footing base)
gamma_concrete = 25  # kN/m^3 (Unit weight of concrete, common value for reinforced)

# Assumed footing dimensions (B=L as per diagram)
B = 10 / 3  # m
L = B  # m

# Eurocode 7 DA1-C1 partial factors
gamma_G = 1.35  # for permanent actions
gamma_Q = 1.50  # for variable actions

# Step 2: Calculate Footing Self-Weight (Characteristic)
Wk_footing = (B * L * D_footing) * gamma_concrete
Gk_total = Gk_structure + Wk_footing

# Step 3: Calculate Design Actions (ULS)
Vd = Gk_total * gamma_G + Qk_v * gamma_Q
Hd = Qk_h * gamma_Q
lever_arm = h_column + D_footing
Md = Hd * lever_arm

# Step 4: Determine Effective Area
# Eccentricity
e = Md / Vd
# Effective dimensions
B_prime = B - 2 * e
A_prime = B_prime * L

# Step 5: Calculate Required ULS Design Bearing Resistance (q_Ed)
q_Ed = Vd / A_prime

# --- Output the results step-by-step ---
print("--- Step-by-step Calculation ---")
print(f"Assumed footing width B = {B:.3f} m")
print(f"1. Total characteristic permanent load Gk_total = Gk_structure + Wk_footing = {Gk_structure} + ({B:.3f} * {L:.3f} * {D_footing} * {gamma_concrete}) = {Gk_total:.2f} kN")
print(f"2. Design vertical load Vd = Gk_total * Y_G + Qk_v * Y_Q = {Gk_total:.2f} * {gamma_G} + {Qk_v} * {gamma_Q} = {Vd:.2f} kN")
print(f"3. Design moment Md = (Qk_h * Y_Q) * (h + D) = ({Qk_h} * {gamma_Q}) * ({h_column} + {D_footing}) = {Md:.2f} kNm")
print(f"4. Eccentricity e = Md / Vd = {Md:.2f} / {Vd:.2f} = {e:.4f} m")
print(f"5. Effective width B' = B - 2*e = {B:.3f} - 2 * {e:.4f} = {B_prime:.4f} m")
print(f"6. Effective area A' = B' * L = {B_prime:.4f} * {L:.3f} = {A_prime:.4f} m^2")
print("\n--- Final Calculation for Required Bearing Resistance (q_Ed) ---")
print(f"q_Ed = Vd / A' = {Vd:.2f} kN / {A_prime:.4f} m^2")
print(f"q_Ed = {q_Ed:.2f} kN/m^2")

# Final answer in the required format
# print(f"\n<<<{q_Ed:.1f}>>>")