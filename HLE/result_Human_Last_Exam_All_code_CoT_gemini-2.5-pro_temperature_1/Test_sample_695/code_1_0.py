import math

def solve_accumulation():
    """
    Calculates the chemical accumulation rate in fish for PFOS and PFOA based on the provided environmental model.
    """

    # --- 1. Define Constants and Parameters ---
    # Environmental parameters
    V = 10000  # L (Volume of freshwater body)
    Qin = 900  # L/d (Inflow rate)
    Qout = 1600 # L/d (Outflow rate)
    foc = 0.1 / 100  # dimensionless (Fraction of organic carbon)
    t = 365  # days (Time)
    C0 = 0   # ng/L (Initial concentration in water, assumed 0 for pristine)

    # Fish parameters (shared)
    Mfish = 1000  # g (Weight of fish)
    Cfood = 100  # ng/g (Concentration of each chemical in food)
    IRfood = 20  # g/day (Ingestion rate of food)
    AFgills = 0.8  # dimensionless (Absorption fraction through gills)
    AFfood = 0.9  # dimensionless (Absorption fraction for food)
    Qgills = 100  # L/day (Gill flow rate)

    # PFOS parameters
    Cin_pfos = 2.6  # ng/L
    half_life_pfos_years = 91  # years
    logKow_pfos = 4.0
    kelim_pfos = 0.069  # days^-1
    Cfish_pfos = 10  # ng/g

    # PFOA parameters
    Cin_pfoa = 211300  # ng/L
    half_life_pfoa_years = 238  # years
    logKow_pfoa = 4.5
    kelim_pfoa = 0.023  # days^-1
    Cfish_pfoa = 10  # ng/g
    
    print("--- Calculation of Chemical Accumulation Rate in Fish ---\n")

    # --- 2. Shared Intermediate Calculations ---
    # Residence rate (r)
    r = Qout / V

    # --- 3. PFOS Calculation ---
    print("--- For PFOS ---")
    # Convert half-life from years to days
    half_life_pfos_days = half_life_pfos_years * 365
    # Degradation rate constant (k)
    k_pfos = math.log(2) / half_life_pfos_days
    # log Koc, Koc, and Kd
    logKoc_pfos = 0.81 * logKow_pfos + 0.01
    Koc_pfos = 10**logKoc_pfos
    Kd_pfos = Koc_pfos * foc
    
    # Water Concentration C(t) for PFOS
    # Denominator of the loading term in the C(t) equation
    denom_ct_pfos = Qout + Qin * (1 + Kd_pfos * foc)
    # The exponential term becomes negligible (close to 0) as t is large
    exp_term_pfos = math.exp(-(k_pfos + r) * t)
    # C(t) equation calculation
    C_water_pfos = (Cin_pfos * Qin / denom_ct_pfos) * (1 - exp_term_pfos) + C0 * exp_term_pfos
    
    print(f"1. Calculate Water Concentration (C_water_pfos) at {t} days:")
    print(f"   log Koc = 0.81 * {logKow_pfos} + 0.01 = {logKoc_pfos:.4f}")
    print(f"   k = ln(2) / ({half_life_pfos_years} * 365) = {k_pfos:.4e} days⁻¹")
    print(f"   r = {Qout} / {V} = {r:.4f} days⁻¹")
    print(f"   C_water_pfos = ({Cin_pfos} * {Qin} / ({Qout} + {Qin} * (1 + {Koc_pfos:.2f} * {foc} * {foc}))) * (1 - e^-(({k_pfos:.4e} + {r:.4f}) * {t}))")
    print(f"   C_water_pfos = {C_water_pfos:.4f} ng/L\n")

    # POP Accumulation for PFOS
    uptake_gills_pfos = C_water_pfos * Qgills * AFgills
    uptake_food_pfos = Cfood * IRfood * AFfood
    elimination_pfos = kelim_pfos * Cfish_pfos * Mfish
    accumulation_pfos = uptake_gills_pfos + uptake_food_pfos - elimination_pfos
    
    print("2. Calculate PFOS Accumulation Rate:")
    print(f"   Accumulation = (C_water × Qgills × AFgills) + (Cfood × IRfood × AFfood) - (kelim × Cfish × Mfish)")
    print(f"   Accumulation_PFOS = ({C_water_pfos:.4f} ng/L * {Qgills} L/d * {AFgills}) + ({Cfood} ng/g * {IRfood} g/d * {AFfood}) - ({kelim_pfos} d⁻¹ * {Cfish_pfos} ng/g * {Mfish} g)")
    print(f"   Accumulation_PFOS = {uptake_gills_pfos:.2f} + {uptake_food_pfos:.2f} - {elimination_pfos:.2f} = {accumulation_pfos:.2f} ng/day\n")

    # --- 4. PFOA Calculation ---
    print("--- For PFOA ---")
    # Convert half-life from years to days
    half_life_pfoa_days = half_life_pfoa_years * 365
    # Degradation rate constant (k)
    k_pfoa = math.log(2) / half_life_pfoa_days
    # log Koc, Koc, and Kd
    logKoc_pfoa = 0.81 * logKow_pfoa + 0.01
    Koc_pfoa = 10**logKoc_pfoa
    Kd_pfoa = Koc_pfoa * foc

    # Water Concentration C(t) for PFOA
    denom_ct_pfoa = Qout + Qin * (1 + Kd_pfoa * foc)
    exp_term_pfoa = math.exp(-(k_pfoa + r) * t)
    C_water_pfoa = (Cin_pfoa * Qin / denom_ct_pfoa) * (1 - exp_term_pfoa) + C0 * exp_term_pfoa
    
    print(f"1. Calculate Water Concentration (C_water_pfoa) at {t} days:")
    print(f"   log Koc = 0.81 * {logKow_pfoa} + 0.01 = {logKoc_pfoa:.4f}")
    print(f"   k = ln(2) / ({half_life_pfoa_years} * 365) = {k_pfoa:.4e} days⁻¹")
    print(f"   r = {Qout} / {V} = {r:.4f} days⁻¹")
    print(f"   C_water_pfoa = ({Cin_pfoa} * {Qin} / ({Qout} + {Qin} * (1 + {Koc_pfoa:.2f} * {foc} * {foc}))) * (1 - e^-(({k_pfoa:.4e} + {r:.4f}) * {t}))")
    print(f"   C_water_pfoa = {C_water_pfoa:.4f} ng/L\n")

    # POP Accumulation for PFOA
    uptake_gills_pfoa = C_water_pfoa * Qgills * AFgills
    uptake_food_pfoa = Cfood * IRfood * AFfood
    elimination_pfoa = kelim_pfoa * Cfish_pfoa * Mfish
    accumulation_pfoa = uptake_gills_pfoa + uptake_food_pfoa - elimination_pfoa
    
    print("2. Calculate PFOA Accumulation Rate:")
    print(f"   Accumulation = (C_water × Qgills × AFgills) + (Cfood × IRfood × AFfood) - (kelim × Cfish × Mfish)")
    print(f"   Accumulation_PFOA = ({C_water_pfoa:.4f} ng/L * {Qgills} L/d * {AFgills}) + ({Cfood} ng/g * {IRfood} g/d * {AFfood}) - ({kelim_pfoa} d⁻¹ * {Cfish_pfoa} ng/g * {Mfish} g)")
    print(f"   Accumulation_PFOA = {uptake_gills_pfoa:.2f} + {uptake_food_pfoa:.2f} - {elimination_pfoa:.2f} = {accumulation_pfoa:.2f} ng/day\n")

    # --- 5. Total Accumulation ---
    total_accumulation = accumulation_pfos + accumulation_pfoa
    print("--- Total Accumulation Rate ---")
    print(f"Total Rate = Accumulation_PFOS + Accumulation_PFOA")
    print(f"Total Rate = {accumulation_pfos:.2f} ng/day + {accumulation_pfoa:.2f} ng/day")
    print(f"Total Chemical Accumulation Rate = {total_accumulation:.2f} ng/day")
    
    return total_accumulation

if __name__ == '__main__':
    final_answer = solve_accumulation()
    # The final answer is requested in a specific format
    print(f"<<<{final_answer:.2f}>>>")
