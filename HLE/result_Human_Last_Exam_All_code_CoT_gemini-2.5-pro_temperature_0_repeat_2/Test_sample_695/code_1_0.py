import math

def calculate_and_print_accumulation(chem_name, params):
    """
    Calculates and prints the water concentration and fish accumulation rate for a given chemical.
    """
    # Extract chemical-specific parameters from the dictionary
    Cin = params['Cin']
    half_life_years = params['half_life']
    log_kow = params['log_kow']
    Cfood = params['Cfood']
    kelim = params['kelim']
    Cfish = params['Cfish']

    # Shared environmental and fish parameters
    V = 10000      # L (Volume of freshwater body)
    Qin = 900      # L/d (Inflow rate)
    Qout = 1600    # L/d (Outflow rate)
    foc = 0.001    # 0.1% (Fraction of organic carbon)
    t = 365        # days
    Qgills = 100   # L/day (Gill flow rate)
    AFgills = 0.8  # (Absorption fraction for gills)
    IRfood = 20    # g/day (Ingestion rate of food)
    AFfood = 0.9   # (Absorption fraction for food)
    Mfish = 1000   # g (Mass of fish)

    # --- Step 1: Calculate intermediate constants ---
    log_koc = 0.81 * log_kow + 0.01
    koc = 10**log_koc
    half_life_days = half_life_years * 365
    k = math.log(2) / half_life_days
    r = Qout / V

    # --- Step 2: Calculate water concentration C(t) at t=365 days ---
    # Using the formula: C(t) = (Cin * Qin / (Qout + Qin * (1 + Koc * foc))) * (1 - e^-(k + r)t)
    # The term with C0 is ignored as C0=0.
    denominator_ct = Qout + Qin * (1 + koc * foc)
    exponent_ct = -(k + r) * t
    C_t = (Cin * Qin / denominator_ct) * (1 - math.exp(exponent_ct))

    # --- Step 3: Calculate POP accumulation rate in fish ---
    # Using the formula: Accumulation = C(t)*Qgills*AFgills + Cfood*IRfood*AFfood - kelim*Cfish*Mfish
    uptake_gills = C_t * Qgills * AFgills
    uptake_food = Cfood * IRfood * AFfood
    elimination = kelim * Cfish * Mfish
    accumulation = uptake_gills + uptake_food - elimination

    # --- Step 4: Print the results with detailed equations ---
    print(f"--- Calculation for {chem_name} ---")
    print(f"Water Concentration C(t) at t={t} days:")
    print(f"C({t}) = ({Cin} * {Qin} / ({Qout} + {Qin} * (1 + {koc:.2f} * {foc}))) * (1 - e^−({k:.4e} + {r:.2f})*{t})")
    print(f"C({t}) = {C_t:.3f} ng/L")
    print("\nPOP Accumulation Rate in Fish:")
    print(f"Accumulation = C(t) * Qgills * AFgills + Cfood * IRfood * AFfood - kelim * Cfish * Mfish")
    print(f"Accumulation = {C_t:.3f} * {Qgills} * {AFgills} + {Cfood} * {IRfood} * {AFfood} - {kelim} * {Cfish} * {Mfish}")
    print(f"Accumulation = {uptake_gills:.2f} + {uptake_food:.2f} - {elimination:.2f} = {accumulation:.2f} ng/day")
    print("-" * (23 + len(chem_name)))
    
    return accumulation

# Define parameters for each chemical in dictionaries
pfos_params = {
    'Cin': 2.6, 'half_life': 91, 'log_kow': 4.0,
    'Cfood': 100, 'kelim': 0.069, 'Cfish': 10
}

pfoa_params = {
    'Cin': 211300, 'half_life': 238, 'log_kow': 4.5,
    'Cfood': 100, 'kelim': 0.023, 'Cfish': 10
}

# Calculate and print results for PFOS
acc_pfos = calculate_and_print_accumulation("PFOS", pfos_params)

print("\n") # Add a newline for separation

# Calculate and print results for PFOA
acc_pfoa = calculate_and_print_accumulation("PFOA", pfoa_params)

final_answer = f"The accumulation rate for PFOS is {acc_pfos:.2f} ng/day and for PFOA is {acc_pfoa:.2f} ng/day."
print(f"\n<<<{final_answer}>>>")