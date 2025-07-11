import math

def calculate_accumulation(chemical_name, Cin, t_half_years, logKow, kelim):
    """
    Calculates the chemical concentration in water and the accumulation rate in fish.
    """
    print(f"--- Calculating for {chemical_name} ---")

    # --- Shared Environmental and Fish Parameters ---
    V_water = 10000  # L
    Qin = 900  # L/d
    Qout = 1600  # L/d
    foc = 0.001  # 0.1% organic carbon
    t = 365  # days
    C0_water = 0 # ng/L, assuming initial concentration is zero

    M_fish = 1000  # g
    C_food = 100  # ng/g
    IR_food = 20  # g/day
    AF_gills = 0.8
    AF_food = 0.9
    Q_gills = 100  # L/day
    C_fish = 10  # ng/g

    # --- Step 1: Calculate intermediate coefficients for water concentration ---

    # Convert half-life from years to days
    t_half_days = t_half_years * 365
    # Calculate decay rate constant (k)
    k = math.log(2) / t_half_days

    # Calculate organic carbon-water partition coefficient (Koc)
    logKoc = 0.81 * logKow + 0.01
    Koc = 10**logKoc

    # Calculate sediment-water partition coefficient (Kd)
    Kd = Koc * foc
    
    # Calculate hydraulic flushing rate (r)
    r = Qout / V_water
    
    # --- Step 2: Calculate the chemical concentration in water C(t) at t=365 days ---
    
    # Using the formula: C(t) = (Cin * Qin / (Qout + Qin * (1 + Kd * foc))) * (1 - e^-(k + r)t) + C0 * e^-(k + r)t
    # Since C0 is 0, the second part of the sum is 0.
    
    css_factor_numerator = Cin * Qin
    css_factor_denominator = Qout + Qin * (1 + Kd * foc)
    css_factor = css_factor_numerator / css_factor_denominator
    
    exponent_val = -(k + r) * t
    
    C_t = css_factor * (1 - math.exp(exponent_val)) + C0_water * math.exp(exponent_val)

    print(f"Calculated water concentration C(365) for {chemical_name}: {C_t:.4f} ng/L")

    # --- Step 3: Calculate the POP accumulation rate in fish ---
    
    # Using the formula: POP accumulation = C(t)*Qgills*AFgills + Cfood*IRfood*AFfood - kelim*Cfish*Mfish
    
    uptake_gills = C_t * Q_gills * AF_gills
    uptake_food = C_food * IR_food * AF_food
    elimination = kelim * C_fish * M_fish
    
    total_accumulation_rate = uptake_gills + uptake_food - elimination

    print("\nFinal Accumulation Equation:")
    print(f"POP Accumulation Rate ({chemical_name}) = C(t) * Qgills * AFgills + Cfood * IRfood * AFfood - kelim * Cfish * Mfish")
    print(f"POP Accumulation Rate ({chemical_name}) = {C_t:.4f} ng/L * {Q_gills} L/day * {AF_gills} + {C_food} ng/g * {IR_food} g/day * {AF_food} - {kelim} days⁻¹ * {C_fish} ng/g * {M_fish} g")
    print(f"POP Accumulation Rate ({chemical_name}) = {uptake_gills:.4f} ng/day + {uptake_food:.4f} ng/day - {elimination:.4f} ng/day")
    print(f"Final Accumulation Rate for {chemical_name}: {total_accumulation_rate:.4f} ng/day\n")
    
    return total_accumulation_rate

# --- Chemical-specific parameters ---
# PFOS
Cin_pfos = 2.6  # ng/L
t_half_pfos = 91  # years
logKow_pfos = 4.0
kelim_pfos = 0.069  # days⁻¹

# PFOA
Cin_pfoa = 211300  # ng/L
t_half_pfoa = 238  # years
logKow_pfoa = 4.5
kelim_pfoa = 0.023  # days⁻¹

# --- Run Calculations ---
pfos_accumulation = calculate_accumulation("PFOS", Cin_pfos, t_half_pfos, logKow_pfos, kelim_pfos)
pfoa_accumulation = calculate_accumulation("PFOA", Cin_pfoa, t_half_pfoa, logKow_pfoa, kelim_pfoa)

pfos_result_str = f'{pfos_accumulation:.4f}'
pfoa_result_str = f'{pfoa_accumulation:.4f}'
# The prompt asks for "the chemical concentration accumulated", but the provided formula gives a rate. We output this calculated rate.
# <<<The accumulation rate for PFOS is 1184.8324 ng/day and for PFOA is 6077131.0264 ng/day>>>