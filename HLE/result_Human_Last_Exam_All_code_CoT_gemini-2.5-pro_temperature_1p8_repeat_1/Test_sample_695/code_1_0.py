import math

def solve_accumulation():
    # --- Step 1: Define all given parameters ---

    # Environmental Parameters
    V = 10000.0  # L (Volume of water body)
    Qin = 900.0   # L/d (Inflow rate)
    Qout = 1600.0 # L/d (Outflow rate)
    foc = 0.001   # 0.1% (Fraction of organic carbon)
    t = 365.0     # days (Time)

    # Fish & Food Parameters
    Mfish = 1000.0   # g (Fish weight)
    IRfood = 20.0    # g/day (Fish ingestion rate)
    Cfood = 100.0    # ng/g (Concentration in food for each chemical)
    AFgills = 0.8    # (Absorption fraction, gills)
    AFfood = 0.9     # (Absorption fraction, food)
    Qgills = 100.0   # L/day (Gill flow rate)
    
    # Assume initial concentration in water C0 is 0
    C0 = 0.0

    # PFOS Parameters
    Cin_pfos = 2.6       # ng/L
    half_life_pfos_y = 91.0  # years
    logKow_pfos = 4.0
    kelim_pfos = 0.069   # days⁻¹
    Cfish_pfos = 10.0    # ng/g

    # PFOA Parameters
    Cin_pfoa = 211300.0  # ng/L
    half_life_pfoa_y = 238.0 # years
    logKow_pfoa = 4.5
    kelim_pfoa = 0.023   # days⁻¹
    Cfish_pfoa = 10.0    # ng/g

    # --- Calculations for PFOS ---
    print("--- Calculating for PFOS ---")

    # Step 2a: Calculate Koc for PFOS
    logKoc_pfos = 0.81 * logKow_pfos + 0.01
    Koc_pfos = 10**logKoc_pfos
    Kd_pfos = Koc_pfos * foc
    
    # Step 2b: Calculate rate constants for PFOS
    half_life_pfos_d = half_life_pfos_y * 365.25 # converting years to days
    k_pfos = math.log(2) / half_life_pfos_d
    r = Qout / V
    k_plus_r_pfos = k_pfos + r

    # Step 2c: Calculate water concentration C(t) for PFOS
    # Denominator from the C(t) equation
    denominator_pfos = Qout + Qin * (1 + Kd_pfos * foc)
    Ct_pfos = (Cin_pfos * Qin / denominator_pfos) * (1 - math.exp(-k_plus_r_pfos * t)) + C0 * math.exp(-k_plus_r_pfos * t)
    print(f"Water concentration of PFOS C(t) at {t} days: {Ct_pfos:.4f} ng/L")

    # Step 3: Calculate POP accumulation rate for PFOS
    uptake_gills_pfos = Ct_pfos * Qgills * AFgills
    uptake_food = Cfood * IRfood * AFfood
    elimination_pfos = kelim_pfos * Cfish_pfos * Mfish
    accumulation_pfos = uptake_gills_pfos + uptake_food - elimination_pfos
    
    print("\nPFOS Accumulation Rate = C(t) * Qgills * AFgills + Cfood * IRfood * AFfood - kelim * Cfish * Mfish")
    print(f"PFOS Accumulation Rate = {Ct_pfos:.4f} * {Qgills} * {AFgills} + {Cfood} * {IRfood} * {AFfood} - {kelim_pfos} * {Cfish_pfos} * {Mfish}")
    print(f"PFOS Accumulation Rate = {uptake_gills_pfos:.2f} + {uptake_food:.2f} - {elimination_pfos:.2f} = {accumulation_pfos:.2f} ng/day")

    # --- Calculations for PFOA ---
    print("\n--- Calculating for PFOA ---")
    
    # Step 2a: Calculate Koc for PFOA
    logKoc_pfoa = 0.81 * logKow_pfoa + 0.01
    Koc_pfoa = 10**logKoc_pfoa
    Kd_pfoa = Koc_pfoa * foc

    # Step 2b: Calculate rate constants for PFOA
    half_life_pfoa_d = half_life_pfoa_y * 365.25 # converting years to days
    k_pfoa = math.log(2) / half_life_pfoa_d
    k_plus_r_pfoa = k_pfoa + r

    # Step 2c: Calculate water concentration C(t) for PFOA
    denominator_pfoa = Qout + Qin * (1 + Kd_pfoa * foc)
    Ct_pfoa = (Cin_pfoa * Qin / denominator_pfoa) * (1 - math.exp(-k_plus_r_pfoa * t)) + C0 * math.exp(-k_plus_r_pfoa * t)
    print(f"Water concentration of PFOA C(t) at {t} days: {Ct_pfoa:.4f} ng/L")

    # Step 3: Calculate POP accumulation rate for PFOA
    uptake_gills_pfoa = Ct_pfoa * Qgills * AFgills
    # uptake_food is the same as for PFOS
    elimination_pfoa = kelim_pfoa * Cfish_pfoa * Mfish
    accumulation_pfoa = uptake_gills_pfoa + uptake_food - elimination_pfoa
    
    print("\nPFOA Accumulation Rate = C(t) * Qgills * AFgills + Cfood * IRfood * AFfood - kelim * Cfish * Mfish")
    print(f"PFOA Accumulation Rate = {Ct_pfoa:.4f} * {Qgills} * {AFgills} + {Cfood} * {IRfood} * {AFfood} - {kelim_pfoa} * {Cfish_pfoa} * {Mfish}")
    print(f"PFOA Accumulation Rate = {uptake_gills_pfoa:.2f} + {uptake_food:.2f} - {elimination_pfoa:.2f} = {accumulation_pfoa:.2f} ng/day")

    # --- Step 4: Sum the results ---
    total_accumulation = accumulation_pfos + accumulation_pfoa
    print("\n--- Total Accumulation ---")
    print(f"Total Accumulation Rate = PFOS Rate + PFOA Rate")
    print(f"Total Accumulation Rate = {accumulation_pfos:.2f} + {accumulation_pfoa:.2f} = {total_accumulation:.2f} ng/day")
    
    return total_accumulation

# Execute the function and store the final answer
final_answer = solve_accumulation()