import math

def calculate_chemical_accumulation():
    """
    Calculates the net accumulation rate of PFOS and PFOA in fish
    based on the provided environmental parameters and equations.
    """
    # --- Step 1: Define Common Environmental and Fish Parameters ---
    V = 10000      # Volume of freshwater body (L)
    Qin = 900      # Inflow rate (L/d)
    Qout = 1600    # Outflow rate (L/d)
    foc = 0.001    # Fraction of organic carbon (0.1%)
    C0 = 0         # Initial concentration in water (ng/L)
    t = 365        # Time (days)

    # Fish parameters
    Mfish = 1000   # Weight of fish (g)
    Cfood = 100    # Concentration in food (ng/g)
    IRfood = 20    # Ingestion rate of food (g/day)
    AFgills = 0.8  # Absorption fraction through gills
    AFfood = 0.9   # Absorption fraction for food
    Qgills = 100   # Gill flow rate (L/day)
    Cfish = 10     # Existing concentration in fish (ng/g)

    # Note on volume: The model assumes constant volume, so V is treated as 10000 L
    # despite Qin != Qout. The dilution rate 'r' is based on outflow.
    r = Qout / V

    # --- Step 2: PFOS Calculation ---
    print("--- PFOS Analysis ---")
    
    # PFOS-specific parameters
    Cin_pfos = 2.6       # Inflow concentration (ng/L)
    t_half_pfos = 91     # Half-life (years)
    log_Kow_pfos = 4.0   # Log octanol-water partition coefficient
    kelim_pfos = 0.069   # Elimination rate constant (days⁻¹)

    # Calculate intermediate values for water concentration of PFOS
    log_Koc_pfos = 0.81 * log_Kow_pfos + 0.01
    Koc_pfos = 10**log_Koc_pfos
    Kd_pfos = Koc_pfos  # Partition coefficient (L/kg)
    k_pfos = math.log(2) / (t_half_pfos * 365.25) # Degradation rate constant (days⁻¹)

    # Calculate PFOS concentration in water at t=365 days using the provided formula
    # C(t) = (Cin * Qin / (Qout + Qin * (1 + Kd * foc))) * (1-e^-(k + r)t) + C0 * e^-(k + r)t
    Ct_pfos_numerator = Cin_pfos * Qin
    Ct_pfos_denominator = Qout + Qin * (1 + Kd_pfos * foc)
    Ct_pfos_ss_part = Ct_pfos_numerator / Ct_pfos_denominator
    exponent_pfos = -(k_pfos + r) * t
    Ct_pfos = Ct_pfos_ss_part * (1 - math.exp(exponent_pfos)) + C0 * math.exp(exponent_pfos)
    
    print(f"Water concentration of PFOS at {t} days is: {Ct_pfos:.4f} ng/L")

    # Calculate net accumulation rate in fish for PFOS
    uptake_gills_pfos = Ct_pfos * Qgills * AFgills
    uptake_food = Cfood * IRfood * AFfood # Same for both chemicals
    elimination_pfos = kelim_pfos * Cfish * Mfish
    accumulation_pfos = uptake_gills_pfos + uptake_food - elimination_pfos
    
    print("\nPFOS net accumulation rate in fish at day 365:")
    print("Equation: C(t) * Qgills * AFgills + Cfood * IRfood * AFfood - kelim * Cfish * Mfish")
    print(f"Calculation: {Ct_pfos:.4f} * {Qgills} * {AFgills} + {Cfood} * {IRfood} * {AFfood} - {kelim_pfos} * {Cfish} * {Mfish}")
    print(f"Result: {uptake_gills_pfos:.2f} + {uptake_food:.2f} - {elimination_pfos:.2f} = {accumulation_pfos:.2f} ng/day")

    print("\n" + "="*50 + "\n")

    # --- Step 3: PFOA Calculation ---
    print("--- PFOA Analysis ---")

    # PFOA-specific parameters
    Cin_pfoa = 211300    # Inflow concentration (ng/L)
    t_half_pfoa = 238    # Half-life (years)
    log_Kow_pfoa = 4.5   # Log octanol-water partition coefficient
    kelim_pfoa = 0.023   # Elimination rate constant (days⁻¹)

    # Calculate intermediate values for water concentration of PFOA
    log_Koc_pfoa = 0.81 * log_Kow_pfoa + 0.01
    Koc_pfoa = 10**log_Koc_pfoa
    Kd_pfoa = Koc_pfoa  # Partition coefficient (L/kg)
    k_pfoa = math.log(2) / (t_half_pfoa * 365.25) # Degradation rate constant (days⁻¹)

    # Calculate PFOA concentration in water at t=365 days
    Ct_pfoa_numerator = Cin_pfoa * Qin
    Ct_pfoa_denominator = Qout + Qin * (1 + Kd_pfoa * foc)
    Ct_pfoa_ss_part = Ct_pfoa_numerator / Ct_pfoa_denominator
    exponent_pfoa = -(k_pfoa + r) * t
    Ct_pfoa = Ct_pfoa_ss_part * (1 - math.exp(exponent_pfoa)) + C0 * math.exp(exponent_pfoa)

    print(f"Water concentration of PFOA at {t} days is: {Ct_pfoa:.4f} ng/L")

    # Calculate net accumulation rate in fish for PFOA
    uptake_gills_pfoa = Ct_pfoa * Qgills * AFgills
    elimination_pfoa = kelim_pfoa * Cfish * Mfish
    accumulation_pfoa = uptake_gills_pfoa + uptake_food - elimination_pfoa

    print("\nPFOA net accumulation rate in fish at day 365:")
    print("Equation: C(t) * Qgills * AFgills + Cfood * IRfood * AFfood - kelim * Cfish * Mfish")
    print(f"Calculation: {Ct_pfoa:.4f} * {Qgills} * {AFgills} + {Cfood} * {IRfood} * {AFfood} - {kelim_pfoa} * {Cfish} * {Mfish}")
    print(f"Result: {uptake_gills_pfoa:.2f} + {uptake_food:.2f} - {elimination_pfoa:.2f} = {accumulation_pfoa:.2f} ng/day")
    
    # Final answer in the requested format
    return accumulation_pfos, accumulation_pfoa

if __name__ == '__main__':
    pfos_result, pfoa_result = calculate_chemical_accumulation()
    print(f"\n<<<{pfos_result:.2f}, {pfoa_result:.2f}>>>")