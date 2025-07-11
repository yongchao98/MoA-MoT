import math

# Step 1: Define Material Properties and Section Geometry
fck = 30.0  # MPa
fyk = 500.0 # MPa
gamma_c = 1.5 # Partial safety factor for concrete
gamma_s = 1.15 # Partial safety factor for steel
Es = 200000.0 # Modulus of elasticity for steel in MPa
eps_cu3 = 0.0035 # Ultimate strain in concrete

# Reinforcement details
d_bar = 20.0 # mm
# Layer 1 (top, compression)
d1 = 50.0 # mm
n1 = 2
# Layer 2 (middle, tension)
d2 = 50.0 + 210.0 # mm
n2 = 2
# Layer 3 (bottom, tension)
d3 = 50.0 + 210.0 + 90.0 # mm
n3 = 3

# Step 2: Calculate Design Strengths and Steel Areas
fcd = fck / gamma_c
fyd = fyk / gamma_s
As_bar = math.pi * (d_bar / 2)**2

As1 = n1 * As_bar # Compression steel
As2 = n2 * As_bar # Tension steel layer 1
As3 = n3 * As_bar # Tension steel layer 2
As_comp = As1
As_tens = As2 + As3

# Step 3: Determine the Neutral Axis Depth (x)
# Equilibrium equation from C = T, assuming all steel yields
# 0.32*fcd*x^2 + 80*fcd*x - fyd*(As_tens - As_comp) = 0
# This is a quadratic equation ax^2 + bx + c = 0
a = 0.32 * fcd
b = 80.0 * fcd
c = -fyd * (As_tens - As_comp)

# Solve for x using the quadratic formula
discriminant = b**2 - 4 * a * c
x = (-b + math.sqrt(discriminant)) / (2 * a)

print("--- Calculation of Neutral Axis (NA) ---")
print(f"Design compressive strength of concrete, fcd = {fcd:.2f} MPa")
print(f"Design yield strength of steel, fyd = {fyd:.2f} MPa")
print(f"Total tension steel area, As,t = {As_tens:.2f} mm^2")
print(f"Compression steel area, As,c = {As_comp:.2f} mm^2")
print(f"From force equilibrium (C=T), we solve for the NA depth 'x'.")
print(f"The equation is: {a:.2f}x^2 + {b:.2f}x + {c:.2f} = 0")
print(f"Solving the quadratic equation gives the neutral axis depth, x = {x:.2f} mm")
print("-" * 20)

# Step 4: Verify Steel Yielding Assumptions
eps_yd = fyd / Es
# Strain in compression steel
eps_sc = eps_cu3 * (x - d1) / x
# Strain in tension steel (layer 2)
eps_st2 = eps_cu3 * (d2 - x) / x
# Strain in tension steel (layer 3)
eps_st3 = eps_cu3 * (d3 - x) / x

print("--- Verification of Steel Strain ---")
print(f"Steel yield strain, ε_yd = {eps_yd:.5f}")
print(f"Strain in compression steel, ε_sc = {eps_sc:.5f}. {'OK, Yields.' if eps_sc >= eps_yd else 'Does not yield.'}")
print(f"Strain in tension steel (layer @{d2}mm), ε_st2 = {eps_st2:.5f}. {'OK, Yields.' if eps_st2 >= eps_yd else 'Does not yield.'}")
print(f"Strain in tension steel (layer @{d3}mm), ε_st3 = {eps_st3:.5f}. {'OK, Yields.' if eps_st3 >= eps_yd else 'Does not yield.'}")
print("-" * 20)

# Step 5: Calculate Moment at Collapse (Mu)
# Depth of rectangular stress block
s = 0.8 * x

# Centroid of the trapezoidal compression area (from top fiber)
zc_num = 50 * s**2 + (s**3 / 3)
zc_den = 100 * s + s**2 / 2
zc = zc_num / zc_den

# Centroid of total tension steel area (from top fiber)
d_tens_centroid = (As2 * d2 + As3 * d3) / (As_tens)

# Forces
Fc = fcd * (100 * s + s**2 / 2)
Fsc = fyd * As_comp

# Lever arms relative to the centroid of tension steel
lever_arm_c = d_tens_centroid - zc
lever_arm_sc = d_tens_centroid - d1

# Moment calculation
Mu = Fc * lever_arm_c + Fsc * lever_arm_sc

# Convert moment to kNm
Mu_kNm = Mu / 1e6

print("--- Moment Calculation ---")
print(f"Depth of rectangular stress block, s = 0.8 * x = {s:.2f} mm")
print(f"Centroid of concrete compression force (from top), zc = {zc:.2f} mm")
print(f"Centroid of tension steel force (from top), d_t = {d_tens_centroid:.2f} mm")
print("\nInternal Forces:")
print(f"Concrete compressive force, Fc = {Fc/1000:.2f} kN")
print(f"Compression steel force, Fsc = {Fsc/1000:.2f} kN")
print("\nLever Arms (relative to tension steel centroid):")
print(f"Lever arm for concrete force, z_c_arm = d_t - zc = {lever_arm_c:.2f} mm")
print(f"Lever arm for compression steel, z_sc_arm = d_t - d1 = {lever_arm_sc:.2f} mm")
print("\nMoment at Collapse (Mu) Equation:")
print(f"Mu = Fc * (d_t - zc) + Fsc * (d_t - d1)")
print(f"Mu = {Fc/1000:.2f} kN * {lever_arm_c:.2f} mm + {Fsc/1000:.2f} kN * {lever_arm_sc:.2f} mm")
print(f"Mu = {Fc * lever_arm_c / 1e6:.2f} kNm + {Fsc * lever_arm_sc / 1e6:.2f} kNm")

print("\n--- Final Answer ---")
print(f"The moment at collapse for the section is: {Mu_kNm:.2f} kNm")

print(f"\n<<<{Mu_kNm:.2f}>>>")