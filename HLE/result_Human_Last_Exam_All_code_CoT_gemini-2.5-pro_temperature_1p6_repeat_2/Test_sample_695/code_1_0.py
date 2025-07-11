import math

def calculate_accumulation():
    """
    Calculates the accumulation rate of PFOS and PFOA in fish based on the provided environmental and chemical data.
    """

    # --- Step 1: Define all given parameters ---

    # Environmental parameters
    V = 10000  # Volume of freshwater body (L)
    Qin = 900    # Inflow rate (L/d)
    Qout = 1600   # Outflow rate (L/d)
    foc = 0.001  # Fraction of organic carbon (0.1%)
    t = 365      # Time (days)
    C0 = 0       # Initial water concentration (assumed, ng/L)

    # Fish parameters
    Mfish = 1000  # Fish weight (g)
    Cfood = 100   # Concentration in food (ng/g) for each chemical
    IRfood = 20   # Ingestion rate of food (g/day)
    AFgills = 0.8 # Absorption fraction through gills
    AFfood = 0.9  # Absorption fraction for food
    Qgills = 100  # Gill flow rate (L/day)
    Cfish_initial = 10 # Given concentration in fish (ng/g)

    # Chemical-specific parameters
    # PFOS
    Cin_pfos = 2.6          # Inflow concentration (ng/L)
    t_half_pfos_years = 91  # Half-life (years)
    logKow_pfos = 4.0       # log Kow
    kelim_pfos = 0.069      # Elimination rate constant (days⁻¹)

    # PFOA
    Cin_pfoa = 211300     # Inflow concentration (ng/L)
    t_half_pfoa_years = 238 # Half-life (years)
    logKow_pfoa = 4.5       # log Kow
    kelim_pfoa = 0.023      # Elimination rate constant (days⁻¹)


    # --- Step 2: Calculation for PFOS ---

    # Calculate intermediate values for C(t) for PFOS
    logKoc_pfos = 0.81 * logKow_pfos + 0.01
    Koc_pfos = 10**logKoc_pfos
    Kd_pfos = Koc_pfos * foc
    k_pfos = math.log(2) / (t_half_pfos_years * 365)
    r = Qout / V

    # Calculate water concentration C(t) for PFOS at t=365 days
    # Using formula: C(t) = (Cin * Qin / (Qout + Qin * (1 + Kd * foc))) * (1 - e^-(k + r)t)
    steady_state_factor_pfos = (Cin_pfos * Qin) / (Qout + Qin * (1 + Kd_pfos * foc))
    exponent_pfos = - (k_pfos + r) * t
    C_water_pfos_365 = steady_state_factor_pfos * (1 - math.exp(exponent_pfos))

    # Calculate POP accumulation rate for PFOS
    # Using formula: Acc = C(t)*Qgills*AFgills + Cfood*IRfood*AFfood - kelim*Cfish*Mfish
    uptake_water_pfos = C_water_pfos_365 * Qgills * AFgills
    uptake_food_pfos = Cfood * IRfood * AFfood
    elimination_pfos = kelim_pfos * Cfish_initial * Mfish
    accumulation_pfos = uptake_water_pfos + uptake_food_pfos - elimination_pfos


    # --- Step 3: Calculation for PFOA ---

    # Calculate intermediate values for C(t) for PFOA
    logKoc_pfoa = 0.81 * logKow_pfoa + 0.01
    Koc_pfoa = 10**logKoc_pfoa
    Kd_pfoa = Koc_pfoa * foc
    k_pfoa = math.log(2) / (t_half_pfoa_years * 365)
    # r is the same for both

    # Calculate water concentration C(t) for PFOA at t=365 days
    steady_state_factor_pfoa = (Cin_pfoa * Qin) / (Qout + Qin * (1 + Kd_pfoa * foc))
    exponent_pfoa = - (k_pfoa + r) * t
    C_water_pfoa_365 = steady_state_factor_pfoa * (1 - math.exp(exponent_pfoa))

    # Calculate POP accumulation rate for PFOA
    uptake_water_pfoa = C_water_pfoa_365 * Qgills * AFgills
    uptake_food_pfoa = Cfood * IRfood * AFfood
    elimination_pfoa = kelim_pfoa * Cfish_initial * Mfish
    accumulation_pfoa = uptake_water_pfoa + uptake_food_pfoa - elimination_pfoa
    
    # --- Step 4: Print results ---

    print("--- PFOS Accumulation Rate in Fish at Day 365 ---")
    print(f"Water concentration of PFOS at 365 days, C(365) = {C_water_pfos_365:.4f} ng/L\n")
    print("Accumulation Rate = (C(365) * Qgills * AFgills) + (Cfood * IRfood * AFfood) - (kelim * Cfish * Mfish)")
    print(f"Accumulation Rate = ({C_water_pfos_365:.4f} ng/L * {Qgills} L/day * {AFgills}) + ({Cfood} ng/g * {IRfood} g/day * {AFfood}) - ({kelim_pfos} days⁻¹ * {Cfish_initial} ng/g * {Mfish} g)")
    print(f"Accumulation Rate = {uptake_water_pfos:.4f} ng/day + {uptake_food_pfos:.4f} ng/day - {elimination_pfos:.4f} ng/day")
    print(f"Net Accumulation Rate (PFOS) = {accumulation_pfos:.4f} ng/day\n")
    
    print("--- PFOA Accumulation Rate in Fish at Day 365 ---")
    print(f"Water concentration of PFOA at 365 days, C(365) = {C_water_pfoa_365:.4f} ng/L\n")
    print("Accumulation Rate = (C(365) * Qgills * AFgills) + (Cfood * IRfood * AFfood) - (kelim * Cfish * Mfish)")
    print(f"Accumulation Rate = ({C_water_pfoa_365:.4f} ng/L * {Qgills} L/day * {AFgills}) + ({Cfood} ng/g * {IRfood} g/day * {AFfood}) - ({kelim_pfoa} days⁻¹ * {Cfish_initial} ng/g * {Mfish} g)")
    print(f"Accumulation Rate = {uptake_water_pfoa:.4f} ng/day + {uptake_food_pfoa:.4f} ng/day - {elimination_pfoa:.4f} ng/day")
    print(f"Net Accumulation Rate (PFOA) = {accumulation_pfoa:.4f} ng/day\n")

    total_accumulation = accumulation_pfos + accumulation_pfoa
    print("--- Total Accumulation Rate ---")
    print(f"Total Net Accumulation Rate (PFOS + PFOA) = {accumulation_pfos:.4f} + {accumulation_pfoa:.4f} = {total_accumulation:.4f} ng/day")
    
    # Return the total value as the final answer
    return total_accumulation

# Execute the function and capture the final answer
final_answer = calculate_accumulation()
# The problem asks for "the chemical concentration accumulated" in plural,
# which has been interpreted as the individual accumulation rates and printed above.
# The single value in the final answer format will be the total accumulation rate.
print(f"\n<<<Total net accumulation rate of {final_answer:.2f} ng/day>>>")
