import math

# --- Given Constants ---

# Environment
V_water = 10000  # L
Qin = 900  # L/d
Qout = 1600  # L/d
foc = 0.001  # 0.1% organic carbon
C0_water = 0 # Assuming pristine environment starts at 0 concentration
t_days = 365 # days

# Fish
M_fish = 1000  # g
C_food = 100  # ng/g of each chemical
IR_food = 20  # g/day
AF_gills = 0.8
AF_food = 0.9
Q_gills = 100  # L/day
C_fish = 10  # ng/g of each chemical at t=365

# Chemical specific constants
# PFOS
Cin_pfos = 2.6  # ng/L
half_life_pfos_years = 91
log_Kow_pfos = 4.0
kelim_pfos = 0.069  # days^-1

# PFOA
Cin_pfoa = 211300  # ng/L
half_life_pfoa_years = 238
log_Kow_pfoa = 4.5
kelim_pfoa = 0.023  # days^-1

def calculate_accumulation(chemical_name, Cin, half_life_years, log_Kow, kelim):
    """
    Calculates the water concentration C(t) and the POP accumulation rate in fish
    for a given chemical.
    """
    print(f"--- Calculating for {chemical_name} ---")

    # --- Step 1: Calculate Water Concentration C(t) at t = 365 days ---

    # Convert half-life from years to days
    t_half_days = half_life_years * 365

    # Calculate first-order decay rate constant (k)
    k = math.log(2) / t_half_days

    # Calculate log Koc and then Koc
    log_Koc = 0.81 * log_Kow + 0.01
    Koc = 10**log_Koc

    # Calculate partition coefficient (Kd)
    Kd = Koc * foc

    # Calculate flushing rate (r)
    r = Qout / V_water

    # Calculate the denominator of the Css term from the provided formula
    # Formula: (Qout + Qin * (1 + Kd * foc))
    denominator_Ct = Qout + Qin * (1 + Kd * foc)

    # Calculate C(t) using the provided formula
    # C(t) = (Cin * Qin / Denominator) * (1 - e^-(k + r)t) + C0 * e^-(k + r)t
    # Since C0 = 0, the second term is 0.
    term1 = (Cin * Qin) / denominator_Ct
    exponent_val = - (k + r) * t_days
    term2 = (1 - math.exp(exponent_val))
    Ct = term1 * term2

    print(f"\n1. Water Concentration Equation for {chemical_name} at {t_days} days:")
    print(f"C({t_days}) = ({Cin:.1f} * {Qin} / ({Qout} + {Qin} * (1 + {Kd:.3f} * {foc}))) * (1 - e^-(({k:.6f} + {r}) * {t_days}))")
    print(f"Calculated Water Concentration C({t_days}): {Ct:.4f} ng/L")

    # --- Step 2: Calculate POP Accumulation Rate in Fish ---
    
    # Uptake from gills (ng/day)
    uptake_gills = Ct * Q_gills * AF_gills
    # Uptake from food (ng/day)
    uptake_food = C_food * IR_food * AF_food
    # Elimination (ng/day)
    elimination = kelim * C_fish * M_fish

    # Net accumulation rate (ng/day)
    accumulation_rate = uptake_gills + uptake_food - elimination

    print(f"\n2. POP Accumulation Rate Equation for {chemical_name}:")
    print(f"Accumulation Rate = C(t) * Qgills * AFgills + Cfood * IRfood * AFfood - kelim * Cfish * Mfish")
    print(f"Accumulation Rate = {Ct:.4f} * {Q_gills} * {AF_gills} + {C_food} * {IR_food} * {AF_food} - {kelim} * {C_fish} * {M_fish}")
    print(f"Calculated Net Accumulation Rate for {chemical_name}: {accumulation_rate:.2f} ng/day\n")

    return accumulation_rate

# --- Perform Calculations ---
accumulation_pfos = calculate_accumulation("PFOS", Cin_pfos, half_life_pfos_years, log_Kow_pfos, kelim_pfos)
accumulation_pfoa = calculate_accumulation("PFOA", Cin_pfoa, half_life_pfoa_years, log_Kow_pfoa, kelim_pfoa)

# Calculate the total accumulation rate for both chemicals
total_accumulation_rate = accumulation_pfos + accumulation_pfoa

print("--- Total Result ---")
print(f"The total net accumulation rate for both chemicals in fish at day 365 is the sum of their individual rates.")
print(f"Total Rate = {accumulation_pfos:.2f} ng/day (PFOS) + {accumulation_pfoa:.2f} ng/day (PFOA)")
print(f"Total Net Accumulation Rate = {total_accumulation_rate:.2f} ng/day")

# Final answer in the required format
# The question asks for "chemical concentration accumulated", which is ambiguous.
# The provided formula calculates a rate (ng/day).
# We provide the total accumulation rate as the final answer.
print(f"\n<<<Total chemical accumulation rate: {total_accumulation_rate:.1f}>>>")