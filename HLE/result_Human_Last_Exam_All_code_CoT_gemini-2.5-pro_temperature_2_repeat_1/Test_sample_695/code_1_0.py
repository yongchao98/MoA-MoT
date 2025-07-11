import math

# --- Given Parameters ---

# Environment
V = 10000      # Volume of freshwater body (L)
Qin = 900      # Inflow rate (L/d)
Qout = 1600    # Outflow rate (L/d)
foc = 0.001    # Organic carbon fraction (0.1%)
C0 = 0         # Initial concentration in water (ng/L) - "pristine"

# Fish
Mfish = 1000   # Average weight of fish (g)
IRfood = 20    # Ingestion rate of food (g/day)
Cfood = 100    # Concentration of each chemical in food (ng/g)
AFgills = 0.8  # Absorption fraction through gills
AFfood = 0.9   # Absorption fraction for food
Qgills = 100   # Gill flow rate (L/day)
Cfish = 10     # Assumed concentration of chemicals in fish for elimination calculation (ng/g)

# Time
t = 365        # days

# --- Chemical-Specific Parameters ---

# PFOS
Cin_pfos = 2.6
t_half_pfos_years = 91
logKow_pfos = 4.0
kelim_pfos = 0.069

# PFOA
Cin_pfoa = 211300
t_half_pfoa_years = 238
logKow_pfoa = 4.5
kelim_pfoa = 0.023


def calculate_accumulation(chemical_name, Cin, t_half_years, logKow, kelim):
    """
    Calculates the net accumulation rate of a chemical in fish.
    
    The function follows the equations and logic provided in the problem description.
    """
    # Step 1: Calculate Koc from log Kow
    # log Koc = 0.81 * log Kow + 0.01
    logKoc = 0.81 * logKow + 0.01
    Koc = 10**logKoc
    
    # Step 2: Calculate terms for the water concentration formula.
    # The term (1 + Kd × foc) is interpreted as (1 + Koc × foc).
    denominator_factor = 1 + Koc * foc
    C_denominator = Qout + Qin * denominator_factor
    
    # Step 3: Calculate environmental degradation rate 'k'
    t_half_days = t_half_years * 365
    k = math.log(2) / t_half_days
    
    # Step 4: Calculate hydraulic removal rate 'r'
    r = Qout / V
    
    # Step 5: Calculate water concentration C(t) at t=365 days
    # C(t) = (Cin * Qin / Denom) * (1 - e^-(k+r)t)
    exponent = -1 * (k + r) * t
    transient_term = 1 - math.exp(exponent) # This will be ~1.0 as the system reaches steady-state
    
    C_t = (Cin * Qin / C_denominator) * transient_term

    # Step 6: Calculate POP accumulation rate in fish
    # POP accumulation = C(t) * Qgills * AFgills + Cfood * IRfood * AFfood - kelim * Cfish * Mfish
    uptake_gills = C_t * Qgills * AFgills
    uptake_food = Cfood * IRfood * AFfood
    elimination = kelim * Cfish * Mfish
    net_accumulation_rate = uptake_gills + uptake_food - elimination
    
    # --- Print Results for this Chemical ---
    print(f"--- {chemical_name} Accumulation Rate Calculation ---")
    print(f"The concentration of {chemical_name} in the water after 365 days, C(365), is {C_t:.4f} ng/L.")
    print("The net accumulation rate is calculated using the equation:")
    print("Accumulation_Rate = C(t) * Qgills * AFgills + Cfood * IRfood * AFfood - kelim * Cfish * Mfish")
    print(f"Accumulation_Rate = {C_t:.4f} * {Qgills} * {AFgills} + {Cfood} * {IRfood} * {AFfood} - {kelim} * {Cfish} * {Mfish}")
    print(f"Resulting in: {net_accumulation_rate:.2f} ng/day = ({uptake_gills:.2f}) + ({uptake_food:.2f}) - ({elimination:.2f})\n")

    return net_accumulation_rate


if __name__ == '__main__':
    acc_rate_pfos = calculate_accumulation(
        chemical_name="PFOS",
        Cin=Cin_pfos,
        t_half_years=t_half_pfos_years,
        logKow=logKow_pfos,
        kelim=kelim_pfos
    )
    
    acc_rate_pfoa = calculate_accumulation(
        chemical_name="PFOA",
        Cin=Cin_pfoa,
        t_half_years=t_half_pfoa_years,
        logKow=logKow_pfoa,
        kelim=kelim_pfoa
    )

    print(f"<<<PFOS accumulation rate: {acc_rate_pfos:.2f} ng/day, PFOA accumulation rate: {acc_rate_pfoa:.2f} ng/day>>>")
