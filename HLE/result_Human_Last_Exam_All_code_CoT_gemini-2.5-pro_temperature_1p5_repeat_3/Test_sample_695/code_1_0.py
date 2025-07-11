import math

def calculate_accumulation():
    # --- Common Parameters ---
    V = 10000  # Volume of freshwater body (L)
    Qin = 900  # Inflow rate (L/d)
    Qout = 1600 # Outflow rate (L/d)
    foc = 0.001  # Fraction of organic carbon (0.1%)
    t = 365  # Time in days
    C0 = 0 # Initial concentration in water (ng/L), assumed pristine

    # --- Fish Parameters ---
    Mfish = 1000  # Weight of fish (g)
    Cfood = 100  # Concentration in food (ng/g)
    IRfood = 20  # Ingestion rate of food (g/day)
    AFgills = 0.8  # Absorption fraction through gills
    AFfood = 0.9  # Absorption fraction for food
    Qgills = 100  # Gill flow rate (L/day)
    Cfish = 10  # Concentration in fish (ng/g) for elimination calculation

    # --- Chemical-Specific Parameters ---
    pfos_params = {
        "name": "PFOS",
        "Cin": 2.6,  # Inflow concentration (ng/L)
        "thalf_years": 91,
        "logKow": 4.0,
        "kelim": 0.069  # Elimination rate constant (days^-1)
    }
    
    pfoa_params = {
        "name": "PFOA",
        "Cin": 211300,  # Inflow concentration (ng/L)
        "thalf_years": 238,
        "logKow": 4.5,
        "kelim": 0.023  # Elimination rate constant (days^-1)
    }

    chemicals = [pfos_params, pfoa_params]
    results = []

    for chem in chemicals:
        # Step 1: Calculate Water Concentration C(t) at t=365 days
        
        # Convert half-life to days and calculate degradation rate k
        thalf_days = chem["thalf_years"] * 365.25
        k = math.log(2) / thalf_days
        
        # Calculate Koc and Kd
        logKoc = 0.81 * chem["logKow"] + 0.01
        Koc = 10**logKoc
        Kd = Koc * foc
        
        # Calculate hydraulic residence rate r
        r = Qout / V
        
        # Calculate C(t) using the provided formula
        # C(t) = (Cin * Qin / (Qout + Qin * (1 + Kd * foc))) * (1 - e^-(k + r)t) + C0 * e^-(k + r)t
        # Since C0 is 0, the second term is 0.
        
        A = (chem["Cin"] * Qin) / (Qout + Qin * (1 + Kd * foc))
        exponent_term = 1 - math.exp(-(k + r) * t)
        C_t = A * exponent_term

        # Step 2: Calculate POP accumulation rate in fish
        
        # POP accumulation = C(t)*Qgills*AFgills + Cfood*IRfood*AFfood - kelim*Cfish*Mfish
        uptake_gills = C_t * Qgills * AFgills
        uptake_food = Cfood * IRfood * AFfood
        elimination = chem["kelim"] * Cfish * Mfish
        
        net_accumulation_rate = uptake_gills + uptake_food - elimination
        results.append(net_accumulation_rate)
        
        print(f"{chem['name']} Accumulation (ng/day):")
        print(f"{net_accumulation_rate:.2f} = {uptake_gills:.2f} (from gills) + {uptake_food:.2f} (from food) - {elimination:.2f} (elimination)")
        print("-" * 20)

    return results

# Run the calculation and store the final answer
final_answers = calculate_accumulation()
print("Final accumulation rates are for PFOS and PFOA respectively.")
print(f"<<<{final_answers[0]:.2f}, {final_answers[1]:.2f}>>>")
