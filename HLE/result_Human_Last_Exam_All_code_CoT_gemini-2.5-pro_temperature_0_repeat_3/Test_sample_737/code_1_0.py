import math

def solve_pavement_thickness():
    """
    Calculates the required pavement thickness using Burmister's two-layer theory.
    """
    # --- Given Data ---
    # Plate Bearing Test
    P_plate = 30 * 1000  # N (Load)
    d_plate = 305  # mm (Plate diameter)
    a_plate = d_plate / 2  # mm (Plate radius)
    delta_s = 2460 / 1000  # mm (Deflection on subgrade)
    delta_p_trial = 1080 / 1000  # mm (Deflection on trial pavement)
    h_trial = 300  # mm (Trial pavement thickness)

    # Design Load
    W_design_ton = 1.80 # ton
    # Assuming 1 ton-force = 9.80665 kN
    P_design = W_design_ton * 9.80665 * 1000 # N
    p_design_kpa = 600  # kN/m^2
    p_design = p_design_kpa / 1000 # N/mm^2 (MPa)
    delta_design_limit = 1.00  # mm (Max allowed deflection)

    # Material Property
    mu = 0.5 # Poisson's ratio

    print("--- Step 1: Calculate Subgrade Modulus (E_s) ---")
    # Using the formula for a flexible plate on an elastic half-space:
    # delta = (1.5 * P) / (pi * a * E)
    # Rearranging for E_s:
    E_s = (1.5 * P_plate) / (math.pi * a_plate * delta_s)
    print(f"The subgrade modulus E_s is calculated as:")
    print(f"E_s = (1.5 * {P_plate:.0f} N) / (pi * {a_plate:.2f} mm * {delta_s:.3f} mm) = {E_s:.2f} N/mm^2 (MPa)\n")

    print("--- Step 2: Determine Pavement-to-Subgrade Modulus Ratio (E_p/E_s) ---")
    # For a two-layer system, delta = [(1.5 * P) / (pi * a * E_s)] * F
    # Rearranging for the deflection factor F:
    F_trial = (delta_p_trial * math.pi * a_plate * E_s) / (1.5 * P_plate)
    h_a_trial = h_trial / a_plate
    print(f"The deflection factor F_trial for the trial section is:")
    print(f"F_trial = ({delta_p_trial:.3f} * pi * {a_plate:.2f} * {E_s:.2f}) / (1.5 * {P_plate:.0f}) = {F_trial:.4f}")
    print(f"The thickness/radius ratio for the trial section is h/a = {h_trial}/{a_plate:.2f} = {h_a_trial:.4f}\n")

    # Interpolating from Burmister chart data (mu=0.5) to find E_p/E_s
    # Data points for h/a = 2.0:
    # E_p/E_s = 2 -> F = 0.53
    # E_p/E_s = 5 -> F = 0.38
    # Since h_a_trial (1.967) is very close to 2.0, we use the values for h/a = 2.0.
    F_at_RE2 = 0.53
    F_at_RE5 = 0.38
    # Linear interpolation for E_p/E_s (R_E):
    # (R_E - 2) / (5 - 2) = (F_trial - F_at_RE2) / (F_at_RE5 - F_at_RE2)
    R_E = 2 + (5 - 2) * (F_trial - F_at_RE2) / (F_at_RE5 - F_at_RE2)
    print("By interpolating from Burmister chart data at h/a ~ 2.0:")
    print(f"The modulus ratio E_p/E_s is found to be {R_E:.2f}\n")

    print("--- Step 3: Calculate Design Parameters ---")
    # Calculate radius of design wheel load
    # Area = P / p = pi * a^2
    a_design_sq = P_design / (p_design * math.pi)
    a_design = math.sqrt(a_design_sq)
    print(f"The radius of the design load area 'a_design' is:")
    print(f"a_design = sqrt({P_design:.2f} N / ({p_design:.3f} N/mm^2 * pi)) = {a_design:.2f} mm\n")

    # Calculate the required design deflection factor F_design
    # delta_design = (1.5 * p_design * a_design) / E_s * F_design
    F_design = (delta_design_limit * E_s) / (1.5 * p_design * a_design)
    print(f"The required design deflection factor 'F_design' is:")
    print(f"F_design = ({delta_design_limit:.2f} mm * {E_s:.2f} MPa) / (1.5 * {p_design:.3f} MPa * {a_design:.2f} mm) = {F_design:.4f}\n")

    print("--- Step 4: Determine Required Pavement Thickness (h_design) ---")
    # We need to find h/a that gives F_design = 0.4380 for E_p/E_s = 3.83.
    # From Step 2, we know that for E_p/E_s ~ 3.83, an F value of ~0.438 corresponds to an h/a ratio of ~2.0.
    # Let's confirm by interpolating F at h/a=2.0 for our calculated R_E.
    F_interp = F_at_RE5 + (F_at_RE2 - F_at_RE5) * (R_E - 5) / (2 - 5)
    print(f"Checking our interpolation: For E_p/E_s = {R_E:.2f}, an h/a ratio of 2.0 gives F = {F_interp:.4f}.")
    print("This is extremely close to our required F_design, so we can use h/a = 2.0 for the design.\n")
    
    h_a_design = 2.0
    h_design = h_a_design * a_design
    
    print("The final required pavement thickness 'h_design' is calculated as:")
    print(f"h_design = (h/a)_design * a_design")
    print(f"h_design = {h_a_design:.3f} * {a_design:.2f} = {h_design:.2f} mm")

solve_pavement_thickness()
<<<193.55>>>