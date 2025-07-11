import math

def calculate_pop_accumulation():
    """
    Calculates the accumulation of PFOS and PFOA in fish based on the provided environmental and biological parameters.
    """

    # --- General and Environmental Parameters ---
    V_water = 10000  # L
    Qin = 900        # L/d
    Qout = 1600      # L/d
    foc = 0.001      # Organic carbon fraction (0.1%)
    t = 365          # days
    C0_water = 0     # Initial water concentration (ng/L), assumed pristine

    # --- Fish Parameters (same for both chemicals) ---
    Mfish = 1000     # g
    Cfood = 100      # ng/g
    IRfood = 20      # g/day
    AFfood = 0.9
    Qgills = 100     # L/day
    AFgills = 0.8
    Cfish = 10       # ng/g

    # --- Chemical-Specific Parameters ---
    chemicals = {
        "PFOS": {
            "Cin": 2.6,         # ng/L
            "t_half_years": 91, # years
            "logKow": 4.0,
            "kelim": 0.069      # days^-1
        },
        "PFOA": {
            "Cin": 211300,      # ng/L
            "t_half_years": 238,# years
            "logKow": 4.5,
            "kelim": 0.023      # days^-1
        }
    }

    # --- Calculations ---
    results = {}
    
    # Calculate flushing rate (same for both)
    r = Qout / V_water

    print("This script calculates the chemical accumulation rate in fish for PFOS and PFOA.")
    print("--------------------------------------------------------------------------------\n")

    for name, params in chemicals.items():
        # Step 1: Calculate water concentration C(t) at t = 365 days
        t_half_days = params["t_half_years"] * 365
        k = math.log(2) / t_half_days
        
        logKoc = 0.81 * params["logKow"] + 0.01
        Koc = 10**logKoc
        Kd = Koc * foc
        
        # Using the provided formula for water concentration C(t)
        # Note: The exponent term (e^-(k+r)t) becomes negligible (approaches 0) due to the large value of (k+r)*t
        # C(t) = (Cin*Qin / (Qout + Qin*(1 + Kd*foc))) * (1-e^-(k+r)t) + C0*e^-(k+r)t
        # For t=365, e^... -> 0, so C(t) simplifies to the steady-state concentration part
        
        denominator = Qout + Qin * (1 + Kd * foc)
        C_ss_term = (params["Cin"] * Qin) / denominator
        exponent_val = -(k + r) * t
        
        # The exponential term is practically zero, but we calculate it for completeness
        exp_term = math.exp(exponent_val)
        C_water_365 = C_ss_term * (1 - exp_term) + C0_water * exp_term

        # Step 2: Calculate the net accumulation rate in fish
        uptake_gills = C_water_365 * Qgills * AFgills
        uptake_food = Cfood * IRfood * AFfood
        elimination = params["kelim"] * Cfish * Mfish
        
        net_accumulation_rate = uptake_gills + uptake_food - elimination
        
        results[name] = net_accumulation_rate

        # --- Print Detailed Steps for Each Chemical ---
        print(f"--- Calculation for {name} ---")
        
        # Print water concentration calculation
        print("\n1. Water Concentration Calculation C(t) at t=365 days:")
        print(f"   logKoc = 0.81 * {params['logKow']} + 0.01 = {logKoc:.4f}")
        print(f"   Koc = 10^{logKoc:.4f} = {Koc:.4f}")
        print(f"   Kd = {Koc:.4f} * {foc} = {Kd:.4f}")
        print(f"   k (degradation rate) = ln(2) / ({params['t_half_years']}*365) = {k:.6f} days⁻¹")
        print(f"   r (flushing rate) = {Qout} / {V_water} = {r} days⁻¹")

        print("\n   Using formula: C(t) = (Cin*Qin / (Qout + Qin*(1 + Kd*foc))) * (1-e^-(k+r)t)")
        print(f"   C(365) = ({params['Cin']} * {Qin} / ({Qout} + {Qin} * (1 + {Kd:.4f} * {foc}))) * (1 - e^-(({k:.6f} + {r}) * {t}))")
        print(f"   Water Concentration of {name} at 365 days: {C_water_365:.4f} ng/L")

        # Print fish accumulation calculation
        print("\n2. Fish Net Accumulation Rate Calculation:")
        print("   Rate = (C_water * Qgills * AFgills) + (C_food * IR_food * AF_food) - (k_elim * C_fish * M_fish)")
        print(f"   Rate = ({C_water_365:.4f} ng/L * {Qgills} L/day * {AFgills}) + ({Cfood} ng/g * {IRfood} g/day * {AFfood}) - ({params['kelim']} days⁻¹ * {Cfish} ng/g * {Mfish} g)")
        print(f"   Rate = ({uptake_gills:.2f} ng/day) + ({uptake_food:.2f} ng/day) - ({elimination:.2f} ng/day)")
        print(f"   >> Final Net Accumulation Rate for {name}: {net_accumulation_rate:.2f} ng/day\n")
        print("--------------------------------------------------------------------------------\n")

    return results

final_results = calculate_pop_accumulation()
pfos_result = final_results.get('PFOS', 0)
pfoa_result = final_results.get('PFOA', 0)

# Final answer in specified format
final_answer_str = f"PFOS accumulation rate: {pfos_result:.2f} ng/day, PFOA accumulation rate: {pfoa_result:.2f} ng/day"
print(f"Final Summary:\n{final_answer_str}")
print(f"<<<{final_answer_str}>>>")