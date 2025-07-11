import math

def solve_accumulation():
    """
    Calculates the concentration of PFOS and PFOA in fish after 365 days.
    """
    # --- Step 0: Define Initial Parameters ---

    # Environment parameters
    V = 10000  # Volume of water body (L)
    Qin = 900  # Inflow rate (L/d)
    Qout = 1600 # Outflow rate (L/d)
    foc = 0.001 # Fraction of organic carbon (0.1%)
    C0 = 0      # Initial concentration in water (ng/L)
    t = 365     # Time (days)

    # Fish parameters
    Mfish = 1000  # Fish weight (g)
    Cfood = 100   # Food concentration (ng/g)
    IRfood = 20   # Food ingestion rate (g/day)
    AFgills = 0.8 # Gill absorption fraction
    AFfood = 0.9  # Food absorption fraction
    Qgills = 100  # Gill flow rate (L/day)
    Cfish_initial = 10 # Initial fish concentration (ng/g)

    # PFOS parameters
    Cin_pfos = 2.6          # Inflow concentration (ng/L)
    half_life_pfos_years = 91
    logKow_pfos = 4.0
    kelim_pfos = 0.069      # Elimination rate (days⁻¹)

    # PFOA parameters
    Cin_pfoa = 211300     # Inflow concentration (ng/L)
    half_life_pfoa_years = 238
    logKow_pfoa = 4.5
    kelim_pfoa = 0.023      # Elimination rate (days⁻¹)

    # --- Step 1: Calculate Water Concentration C(t) for PFOS and PFOA ---
    print("--- Part 1: Calculating Water Concentration after 365 days ---\n")

    # Common water parameters
    r = Qout / V  # Residence rate (d⁻¹)

    # PFOS Water Calculation
    half_life_pfos_days = half_life_pfos_years * 365
    k_pfos = math.log(2) / half_life_pfos_days
    logKoc_pfos = 0.81 * logKow_pfos + 0.01
    Koc_pfos = 10**logKoc_pfos
    
    # Using the provided formula for C(t)
    water_denom_pfos = Qout + Qin * (1 + Koc_pfos * foc)
    water_factor_pfos = (Cin_pfos * Qin) / water_denom_pfos
    k_plus_r_pfos = k_pfos + r
    time_factor_pfos = (1 - math.exp(-k_plus_r_pfos * t))
    C_water_pfos = water_factor_pfos * time_factor_pfos + C0 * math.exp(-k_plus_r_pfos * t)
    
    print(f"PFOS water concentration (ng/L) = ({Cin_pfos} * {Qin} / ({Qout} + {Qin} * (1 + {Koc_pfos:.2f} * {foc}))) * (1 - e^(-({k_pfos:.6f} + {r}) * {t}))")
    print(f"Calculated PFOS concentration in water C(365): {C_water_pfos:.4f} ng/L\n")

    # PFOA Water Calculation
    half_life_pfoa_days = half_life_pfoa_years * 365
    k_pfoa = math.log(2) / half_life_pfoa_days
    logKoc_pfoa = 0.81 * logKow_pfoa + 0.01
    Koc_pfoa = 10**logKoc_pfoa
    
    water_denom_pfoa = Qout + Qin * (1 + Koc_pfoa * foc)
    water_factor_pfoa = (Cin_pfoa * Qin) / water_denom_pfoa
    k_plus_r_pfoa = k_pfoa + r
    time_factor_pfoa = (1 - math.exp(-k_plus_r_pfoa * t))
    C_water_pfoa = water_factor_pfoa * time_factor_pfoa + C0 * math.exp(-k_plus_r_pfoa * t)

    print(f"PFOA water concentration (ng/L) = ({Cin_pfoa} * {Qin} / ({Qout} + {Qin} * (1 + {Koc_pfoa:.2f} * {foc}))) * (1 - e^(-({k_pfoa:.6f} + {r}) * {t}))")
    print(f"Calculated PFOA concentration in water C(365): {C_water_pfoa:.4f} ng/L\n")
    
    # --- Step 2: Calculate Fish Concentration C_fish(t) for PFOS and PFOA ---
    print("--- Part 2: Calculating Fish Concentration after 365 days ---\n")
    
    # PFOS Fish Calculation
    Uptake_rate_pfos = (C_water_pfos * Qgills * AFgills) + (Cfood * IRfood * AFfood)
    A_pfos = Uptake_rate_pfos / Mfish
    B_pfos = kelim_pfos
    C_ss_pfos = A_pfos / B_pfos
    exp_term_pfos = math.exp(-B_pfos * t)
    C_fish_pfos = C_ss_pfos * (1 - exp_term_pfos) + Cfish_initial * exp_term_pfos

    print("Final PFOS concentration in fish is calculated as: C_ss * (1 - e^(-kelim * t)) + C_initial * e^(-kelim * t)")
    print(f"Final PFOS concentration (ng/g) = {C_ss_pfos:.4f} * (1 - e^(-{B_pfos} * {t})) + {Cfish_initial} * e^(-{B_pfos} * {t})")
    print(f"Resulting PFOS concentration in fish: {C_fish_pfos:.4f} ng/g\n")

    # PFOA Fish Calculation
    Uptake_rate_pfoa = (C_water_pfoa * Qgills * AFgills) + (Cfood * IRfood * AFfood)
    A_pfoa = Uptake_rate_pfoa / Mfish
    B_pfoa = kelim_pfoa
    C_ss_pfoa = A_pfoa / B_pfoa
    exp_term_pfoa = math.exp(-B_pfoa * t)
    C_fish_pfoa = C_ss_pfoa * (1 - exp_term_pfoa) + Cfish_initial * exp_term_pfoa

    print("Final PFOA concentration in fish is calculated as: C_ss * (1 - e^(-kelim * t)) + C_initial * e^(-kelim * t)")
    print(f"Final PFOA concentration (ng/g) = {C_ss_pfoa:.4f} * (1 - e^(-{B_pfoa} * {t})) + {Cfish_initial} * e^(-{B_pfoa} * {t})")
    print(f"Resulting PFOA concentration in fish: {C_fish_pfoa:.4f} ng/g")

    return C_fish_pfos, C_fish_pfoa

if __name__ == '__main__':
    pfos_conc, pfoa_conc = solve_accumulation()
    # The final answer format required by the user prompt
    print(f"<<<{pfos_conc:.2f}, {pfoa_conc:.2f}>>>")
