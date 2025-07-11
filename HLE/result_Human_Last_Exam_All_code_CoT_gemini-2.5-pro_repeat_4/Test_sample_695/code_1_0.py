import math

def calculate_final_fish_concentration():
    """
    This function calculates the final concentration of PFOS and PFOA in fish after 365 days
    based on a differential equation model for bioaccumulation. It interprets the problem
    by first finding the steady-state water concentration, then using that to model the
    change in fish concentration over the specified period.
    """

    # --- Step 1: Define Constants ---
    # Environment
    Qin = 900.0      # L/d (Inflow rate)
    Qout = 1600.0    # L/d (Outflow rate)
    foc = 0.1 / 100.0 # dimensionless fraction of organic carbon
    
    # Fish
    Mfish = 1000.0   # g (Fish weight)
    IRfood = 20.0    # g/day (Ingestion rate of food)
    Cfood = 100.0    # ng/g (Concentration in food)
    AFgills = 0.8    # Absorption fraction from gills
    AFfood = 0.9     # Absorption fraction from food
    Qgills = 100.0   # L/day (Gill flow rate)
    Cfish_initial = 10.0 # ng/g, initial concentration in fish C_fish(0)
    
    # Time
    t = 365.0        # days

    # Chemical properties dictionaries
    pfos_params = {
        'name': 'PFOS',
        'Cin': 2.6,      # ng/L
        'logKow': 4.0,
        'kelim': 0.069   # days⁻¹
    }
    pfoa_params = {
        'name': 'PFOA',
        'Cin': 211300.0, # ng/L
        'logKow': 4.5,
        'kelim': 0.023   # days⁻¹
    }
    chemicals = [pfos_params, pfoa_params]

    final_concentrations = {}

    # --- Step 2: Loop through each chemical ---
    for chem in chemicals:
        print(f"\n--- Calculating for {chem['name']} ---")

        # a. Calculate Koc from logKow using log Koc = 0.81 * log Kow + 0.01
        logKoc = 0.81 * chem['logKow'] + 0.01
        Koc = 10**logKoc
        
        # b. Calculate steady-state water concentration (C_water_ss)
        # Using the steady-state part of the provided C(t) equation: Css = (Cin*Qin)/(Qout + Qin*(1+Kd*foc))
        # We assume Kd is numerically equal to the calculated Koc.
        C_water_ss_numerator = chem['Cin'] * Qin
        C_water_ss_denominator = Qout + Qin * (1 + Koc * foc)
        C_water_ss = C_water_ss_numerator / C_water_ss_denominator
        print(f"Calculated steady-state water concentration (C_water): {C_water_ss:.4f} ng/L")

        # c. Calculate total constant uptake rate (in ng/day) from gills and food
        uptake_gills = C_water_ss * Qgills * AFgills
        uptake_food = Cfood * IRfood * AFfood
        total_uptake_rate = uptake_gills + uptake_food
        
        # d. Solve the bioaccumulation differential equation for Cfish(t)
        # Cfish(t) = (SteadyStateFishConc) * (1 - e^(-kelim*t)) + Cfish(0) * e^(-kelim*t)
        # where SteadyStateFishConc = (TotalUptakeRate / FishMass) / EliminationConstant
        kelim = chem['kelim']
        
        # Components of the solution
        steady_state_fish_conc = (total_uptake_rate / Mfish) / kelim
        exp_term = math.exp(-kelim * t)
        
        Cfish_final = steady_state_fish_conc * (1 - exp_term) + Cfish_initial * exp_term
        
        final_concentrations[chem['name']] = Cfish_final
        
        # Print the final equation with numbers as requested
        print("\nCalculating final fish concentration (C_fish) at t=365 days:")
        print("Equation: C_fish(t) = C_fish_ss * (1 - e^(-k_elim * t)) + C_fish(0) * e^(-k_elim * t)")
        print(f"where C_fish_ss (steady-state fish conc.) = (Uptake / Mass) / k_elim = ({total_uptake_rate:.2f} / {Mfish}) / {kelim} = {steady_state_fish_conc:.2f} ng/g")
        print(f"C_fish(365) = {steady_state_fish_conc:.2f} * (1 - e^(-{kelim} * {t})) + {Cfish_initial} * e^(-{kelim} * {t})")
        print(f"Final Concentration for {chem['name']} after {t} days: {Cfish_final:.2f} ng/g")

    # --- Step 3: Present the final answer ---
    pfos_final_conc = final_concentrations['PFOS']
    pfoa_final_conc = final_concentrations['PFOA']
    
    print("\n" + "="*50)
    print("\nSUMMARY OF RESULTS:")
    print(f"The final concentration of PFOS in fish after 365 days is {pfos_final_conc:.2f} ng/g.")
    print(f"The final concentration of PFOA in fish after 365 days is {pfoa_final_conc:.2f} ng/g.")


# Execute the calculation
calculate_final_fish_concentration()