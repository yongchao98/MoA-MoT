import math

def calculate_accumulation():
    """
    Calculates the water concentration and subsequent net accumulation rate in fish
    for PFOS and PFOA based on the provided environmental parameters and equations.
    """
    # --- Environment and Fish Parameters (shared for both chemicals) ---
    vol_water = 10000  # L
    inflow_rate_Qin = 900  # L/d
    outflow_rate_Qout = 1600  # L/d
    frac_organic_carbon_foc = 0.001  # 0.1%
    time_t = 365  # days
    fish_mass_Mfish = 1000  # g
    food_conc_Cfood = 100  # ng/g
    ingestion_rate_IRfood = 20  # g/day
    absorption_frac_gills_AFgills = 0.8
    absorption_frac_food_AFfood = 0.9
    gill_flow_rate_Qgills = 100  # L/day
    initial_fish_conc_Cfish = 10  # ng/g
    initial_water_conc_C0 = 0 # ng/L, assumed for pristine environment

    # --- Chemical-Specific Parameters ---
    # PFOS
    cin_pfos = 2.6  # ng/L
    half_life_years_pfos = 91  # years
    log_kow_pfos = 4.0
    kelim_pfos = 0.069  # days^-1

    # PFOA
    cin_pfoa = 211300  # ng/L
    half_life_years_pfoa = 238  # years
    log_kow_pfoa = 4.5
    kelim_pfoa = 0.023  # days^-1

    chemicals = {
        "PFOS": {
            "Cin": cin_pfos,
            "t_half_years": half_life_years_pfos,
            "log_Kow": log_kow_pfos,
            "kelim": kelim_pfos
        },
        "PFOA": {
            "Cin": cin_pfoa,
            "t_half_years": half_life_years_pfoa,
            "log_Kow": log_kow_pfoa,
            "kelim": kelim_pfoa
        }
    }

    results = {}

    for name, params in chemicals.items():
        # --- Step 1: Calculate Water Concentration C(t) ---
        log_koc = 0.81 * params["log_Kow"] + 0.01
        koc = 10**log_koc
        kd = koc * frac_organic_carbon_foc
        
        half_life_days = params["t_half_years"] * 365
        k_degradation = math.log(2) / half_life_days
        r_hydraulic = outflow_rate_Qout / vol_water
        
        # Using the provided water concentration equation
        steady_state_numerator = params["Cin"] * inflow_rate_Qin
        steady_state_denominator = outflow_rate_Qout + inflow_rate_Qin * (1 + kd * frac_organic_carbon_foc)
        steady_state_conc = steady_state_numerator / steady_state_denominator
        
        rate_term = k_degradation + r_hydraulic
        exp_term = math.exp(-rate_term * time_t)
        
        c_water_t = steady_state_conc * (1 - exp_term) + initial_water_conc_C0 * exp_term
        
        # --- Step 2: Calculate POP Accumulation Rate ---
        uptake_gills = c_water_t * gill_flow_rate_Qgills * absorption_frac_gills_AFgills
        uptake_food = food_conc_Cfood * ingestion_rate_IRfood * absorption_frac_food_AFfood
        elimination = params["kelim"] * initial_fish_conc_Cfish * fish_mass_Mfish
        
        net_accumulation = uptake_gills + uptake_food - elimination
        results[name] = {
            "c_water_t": c_water_t,
            "net_accumulation": net_accumulation
        }

    # --- Step 3: Print the Final Results with Equations ---
    pfos_c_water = results["PFOS"]["c_water_t"]
    pfos_accumulation = results["PFOS"]["net_accumulation"]
    
    pfoa_c_water = results["PFOA"]["c_water_t"]
    pfoa_accumulation = results["PFOA"]["net_accumulation"]

    print("--- PFOS Accumulation Calculation ---")
    print(f"The net accumulation rate of PFOS in fish after {time_t} days is calculated as follows:")
    print(f"Accumulation = (C(t) * Qgills * AFgills) + (Cfood * IRfood * AFfood) - (kelim * Cfish * Mfish)")
    print(f"Accumulation = ({pfos_c_water:.3f} ng/L * {gill_flow_rate_Qgills} L/day * {absorption_frac_gills_AFgills}) + ({food_conc_Cfood} ng/g * {ingestion_rate_IRfood} g/day * {absorption_frac_food_AFfood}) - ({kelim_pfos} days⁻¹ * {initial_fish_conc_Cfish} ng/g * {fish_mass_Mfish} g)")
    print(f"Result: {pfos_accumulation:.2f} ng/day\n")

    print("--- PFOA Accumulation Calculation ---")
    print(f"The net accumulation rate of PFOA in fish after {time_t} days is calculated as follows:")
    print(f"Accumulation = (C(t) * Qgills * AFgills) + (Cfood * IRfood * AFfood) - (kelim * Cfish * Mfish)")
    print(f"Accumulation = ({pfoa_c_water:.3f} ng/L * {gill_flow_rate_Qgills} L/day * {absorption_frac_gills_AFgills}) + ({food_conc_Cfood} ng/g * {ingestion_rate_IRfood} g/day * {absorption_frac_food_AFfood}) - ({kelim_pfoa} days⁻¹ * {initial_fish_conc_Cfish} ng/g * {fish_mass_Mfish} g)")
    print(f"Result: {pfoa_accumulation:.2f} ng/day")
    
    # Returning final values for potential extraction
    final_answer = {
        "PFOS_accumulation_ng_per_day": pfos_accumulation,
        "PFOA_accumulation_ng_per_day": pfoa_accumulation
    }
    # This print is for machine-readability if needed, based on prompt format.
    # We will use the more human-readable printouts above as the main answer.
    # print(f'<<<{final_answer}>>>')

# Execute the calculation and print the results
calculate_accumulation()