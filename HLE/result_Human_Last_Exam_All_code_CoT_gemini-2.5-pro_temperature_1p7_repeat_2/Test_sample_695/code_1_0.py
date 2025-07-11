import math

def calculate_final_fish_concentration():
    """
    Calculates the final steady-state concentration of PFOS and PFOA in fish
    based on the provided environmental and biological parameters and equations.
    """

    # --- Given Common Parameters ---
    # Environment
    V = 10000.0  # L
    Qin = 900.0  # L/d
    Qout = 1600.0 # L/d
    foc = 0.001 # 0.1% as a fraction
    t = 365.0 # days

    # Fish
    Mfish = 1000.0 # g
    Cfood = 100.0 # ng/g
    IRfood = 20.0 # g/day
    AFgills = 0.8
    AFfood = 0.9
    Qgills = 100.0 # L/day

    # --- Chemical-Specific Parameters ---
    chemicals = {
        "PFOS": {
            "Cin": 2.6, # ng/L
            "half_life_yr": 91.0,
            "logKow": 4.0,
            "kelim": 0.069 # days^-1
        },
        "PFOA": {
            "Cin": 211300.0, # ng/L
            "half_life_yr": 238.0,
            "logKow": 4.5,
            "kelim": 0.023 # days^-1
        }
    }

    print("This script calculates the final concentration of PFOS and PFOA in fish.\n")
    print("Step 1: Calculate the steady-state chemical concentration in the water, C_water.")
    print("Step 2: Calculate the total uptake rate (gills + food) for the fish.")
    print("Step 3: Assuming steady state (Uptake = Elimination), solve for the final fish concentration, C_fish.\n")
    
    final_concentrations = {}

    for name, params in chemicals.items():
        print(f"--- Calculations for {name} ---")

        # --- Step 1: Calculate C_water ---
        logKoc = 0.81 * params["logKow"] + 0.01
        Koc = 10**logKoc
        # Per the flawed formula C(t), we assume Kd=Koc for the calculation.
        Kd = Koc
        
        # Steady-state water concentration C_water = (Cin * Qin) / (Qout + Qin * (1 + Kd * foc))
        c_water_numerator = params["Cin"] * Qin
        c_water_denominator = Qout + Qin * (1 + Kd * foc)
        C_water = c_water_numerator / c_water_denominator

        print(f"1. Water Concentration (C_water):")
        print(f"   logKoc = 0.81 * {params['logKow']} + 0.01 = {logKoc:.3f}")
        print(f"   Koc = 10^{logKoc:.3f} = {Koc:.2f} L/kg")
        print(f"   C_water = ({params['Cin']} ng/L * {Qin} L/d) / ({Qout} L/d + {Qin} L/d * (1 + {Koc:.2f} L/kg * {foc}))")
        print(f"   C_water = {c_water_numerator:.2f} / {c_water_denominator:.2f} = {C_water:.4f} ng/L\n")

        # --- Step 2: Calculate Total Uptake Rate ---
        uptake_gills = C_water * Qgills * AFgills
        uptake_food = Cfood * IRfood * AFfood
        total_uptake_ng_per_day = uptake_gills + uptake_food
        
        print("2. Total Uptake Rate:")
        print(f"   Uptake (gills) = {C_water:.4f} ng/L * {Qgills} L/d * {AFgills} = {uptake_gills:.2f} ng/day")
        print(f"   Uptake (food) = {Cfood} ng/g * {IRfood} g/d * {AFfood} = {uptake_food:.2f} ng/day")
        print(f"   Total Uptake = {uptake_gills:.2f} + {uptake_food:.2f} = {total_uptake_ng_per_day:.2f} ng/day\n")

        # --- Step 3: Calculate Final Fish Concentration ---
        # At steady state, Uptake = Elimination, where Elimination = kelim * C_fish * M_fish
        # Therefore, C_fish = Total_Uptake / (kelim * M_fish)
        C_fish_final = total_uptake_ng_per_day / (params["kelim"] * Mfish)
        
        print("3. Final Fish Concentration (C_fish):")
        print(f"   At steady-state, Uptake Rate = Elimination Rate")
        print(f"   {total_uptake_ng_per_day:.2f} ng/day = {params['kelim']} days⁻¹ * C_fish * {Mfish} g")
        print(f"   C_fish = {total_uptake_ng_per_day:.2f} / ({params['kelim']} * {Mfish})")
        print(f"   Final C_fish ({name}) = {C_fish_final:.4f} ng/g\n")
        
        final_concentrations[name] = C_fish_final

    return final_concentrations

# Run the calculation and print the final answer
final_results = calculate_final_fish_concentration()

print("--- Summary of Final Accumulated Concentrations ---")
print(f"The final accumulated concentration for PFOS is {final_results['PFOS']:.2f} ng/g.")
print(f"The final accumulated concentration for PFOA is {final_results['PFOA']:.2f} ng/g.")
print("<<<")
print(f"PFOS: {final_results['PFOS']:.2f} ng/g")
print(f"PFOA: {final_results['PFOA']:.2f} ng/g")
print(">>>")
