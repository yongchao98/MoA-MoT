import math

def calculate_fish_concentration(chem_name, Cin, half_life_years, log_Kow, kelim, Cfish_initial):
    """
    Calculates the final concentration of a chemical in fish after 365 days.
    """
    # Step 1: Define constants and state assumptions
    # Environmental parameters
    V = 10000  # L
    Qin = 900  # L/d
    Qout = 1600  # L/d
    foc = 0.001  # 0.1% organic carbon
    C0 = 0 # ng/L, pristine environment
    
    # Fish parameters
    Mfish = 1000  # g
    Cfood = 100  # ng/g
    IRfood = 20  # g/day
    AFgills = 0.8
    AFfood = 0.9
    Qgills = 100  # L/day

    # Time
    t = 365  # days
    
    # Assumption: The water volume V is assumed to be constant.
    # The flushing rate 'r' is determined by the outflow rate.
    r = Qout / V # d^-1

    # Step 2: Calculate intermediate values for water concentration
    half_life_days = half_life_years * 365
    k = math.log(2) / half_life_days  # d^-1

    log_Koc = 0.81 * log_Kow + 0.01
    Koc = 10**log_Koc
    Kd = Koc * foc
    
    # Step 3: Calculate the water concentration C(t) at t=365 days
    # Using the formula: C(t) = (Cin * Qin / (Qout + Qin * (1 + Kd * foc))) * (1 - e^-(k + r)t) + C0 * e^-(k + r)t
    # Since t=365 and r is large, the exp(-(k+r)t) term approaches 0, so C(t) approaches steady state.
    water_conc_ss_numerator = Cin * Qin
    water_conc_ss_denominator = Qout + Qin * (1 + Kd * foc)
    C_water_ss = water_conc_ss_numerator / water_conc_ss_denominator
    
    # Calculate C(t) explicitly for accuracy
    exponent_term_water = math.exp(-(k + r) * t)
    C_water_t = C_water_ss * (1 - exponent_term_water) + C0 * exponent_term_water
    
    # Step 4: Calculate the final fish concentration
    # The accumulation formula represents d(Amount)/dt. We solve the derived ODE for C_fish(t).
    # dC_fish/dt = (Uptake / Mfish) - kelim * C_fish
    # Solution: C_fish(t) = (Css_fish) * (1 - e^(-kelim*t)) + Cfish_initial * e^(-kelim*t)
    # where Css_fish = (Uptake / Mfish) / kelim
    
    # Calculate total uptake rate
    uptake_gills = C_water_t * Qgills * AFgills
    uptake_food = Cfood * IRfood * AFfood
    total_uptake = uptake_gills + uptake_food
    
    # Calculate the steady-state fish concentration that would be reached
    Cfish_ss = (total_uptake / Mfish) / kelim
    
    # Calculate the concentration at t=365
    exponent_term_fish = math.exp(-kelim * t)
    Cfish_final = Cfish_ss * (1 - exponent_term_fish) + Cfish_initial * exponent_term_fish
    
    # Step 5: Print the results as requested
    print(f"--- Calculation for {chem_name} ---")
    print("This calculation assumes a constant water volume of 10000 L.")
    print("\n1. Final Water Concentration Calculation:")
    print(f"C_water({t} days) = ({Cin} ng/L * {Qin} L/d / ({Qout} L/d + {Qin} L/d * (1 + {Kd:.2f} * {foc}))) * (1 - e^-(({k:.6f} + {r}) d⁻¹ * {t} days))")
    print(f"Final Water Concentration = {C_water_t:.4f} ng/L")

    print("\n2. Final Fish Concentration Calculation:")
    print(f"Total Uptake Rate = (Water Conc. * Gill Flow * Gill Abs.) + (Food Conc. * Ingestion Rate * Food Abs.)")
    print(f"Total Uptake Rate = ({C_water_t:.4f} ng/L * {Qgills} L/day * {AFgills}) + ({Cfood} ng/g * {IRfood} g/day * {AFfood}) = {total_uptake:.2f} ng/day")
    print("\nC_fish(t) = ((Total_Uptake / M_fish) / k_elim) * (1 - e^(-k_elim * t)) + C_fish_initial * e^(-k_elim * t)")
    print(f"C_fish({t}) = (({total_uptake:.2f} ng/day / {Mfish} g) / {kelim} d⁻¹) * (1 - e^(-{kelim} d⁻¹ * {t} days)) + {Cfish_initial} ng/g * e^(-{kelim} d⁻¹ * {t} days)")
    print(f"\nFinal {chem_name} concentration in fish: {Cfish_final:.4f} ng/g\n")
    return Cfish_final

# --- Main execution ---
# Parameters for PFOS
pfos_params = {
    "chem_name": "PFOS",
    "Cin": 2.6, # ng/L
    "half_life_years": 91,
    "log_Kow": 4.0,
    "kelim": 0.069, # days^-1
    "Cfish_initial": 10 # ng/g
}

# Parameters for PFOA
pfoa_params = {
    "chem_name": "PFOA",
    "Cin": 211300, # ng/L
    "half_life_years": 238,
    "log_Kow": 4.5,
    "kelim": 0.023, # days^-1
    "Cfish_initial": 10 # ng/g
}

final_pfos_conc = calculate_fish_concentration(**pfos_params)
final_pfoa_conc = calculate_fish_concentration(**pfoa_params)

# Final answer in specified format
# The question asks for "the chemical concentration", which is plural, referring to both chemicals.
# We present both calculated concentrations as the answer.
print(f"<<<Final PFOS Concentration: {final_pfos_conc:.4f} ng/g, Final PFOA Concentration: {final_pfoa_conc:.4f} ng/g>>>")
