import math

# Step 1: Define all given parameters

# Environmental Parameters
V = 10000  # Volume of freshwater body (L)
Qin = 900  # Inflow rate (L/d)
Qout = 1600  # Outflow rate (L/d)
foc = 0.001  # Organic carbon fraction (0.1%)
C0 = 0  # Initial water concentration (ng/L), assumed for "pristine" environment
t = 365  # Time (days)

# Fish Parameters
Mfish = 1000  # Fish weight (g)
Cfood = 100  # Food concentration (ng/g of each chemical)
IRfood = 20  # Ingestion rate (g/day)
AFgills = 0.8  # Gill absorption fraction
AFfood = 0.9  # Food absorption fraction
Qgills = 100  # Gill flow rate (L/day)
Cfish = 10  # Fish concentration (ng/g of each chemical)

# PFOS Specific Parameters
Cin_pfos = 2.6  # Inflow concentration of PFOS (ng/L)
half_life_pfos_years = 91  # Half-life of PFOS (years)
log_Kow_pfos = 4.0  # Log Kow of PFOS
kelim_pfos = 0.069  # Elimination rate constant for PFOS (days⁻¹)

# PFOA Specific Parameters
Cin_pfoa = 211300  # Inflow concentration of PFOA (ng/L)
half_life_pfoa_years = 238  # Half-life of PFOA (years)
log_Kow_pfoa = 4.5  # Log Kow of PFOA
kelim_pfoa = 0.023  # Elimination rate constant for PFOA (days⁻¹)

def calculate_water_concentration(Cin, half_life_years, log_Kow):
    """
    Calculates the chemical concentration in water at t=365 days.
    Uses the equation: C(t) = (Cin * Qin / (Qout + Qin * (1 + Kd * foc))) * (1-e^-(k + r)t)
    """
    # Step 2a: Convert half-life to days and calculate degradation rate (k)
    half_life_days = half_life_years * 365.25
    k = math.log(2) / half_life_days

    # Step 2b: Calculate Koc and Kd
    log_Koc = 0.81 * log_Kow + 0.01
    Koc = 10**log_Koc
    Kd = Koc * foc
    
    # Step 2c: Calculate hydraulic flushing rate (r)
    r = Qout / V
    
    # Step 2d: Calculate C(t)
    # The term (k + r) * t will be large, making e^-(k+r)t very close to 0
    # and (1 - e^-(k+r)t) very close to 1.
    exponent_term = math.exp(-(k + r) * t)
    
    # Per user's formula, the term is (Qout + Qin * (1 + Kd * foc))
    denominator = Qout + Qin * (1 + Kd * foc)
    
    concentration_t = (Cin * Qin / denominator) * (1 - exponent_term) + C0 * exponent_term
    return concentration_t

def calculate_accumulation_rate(C_t, kelim):
    """
    Calculates the POP accumulation rate in fish at t=365 days.
    Uses the equation: rate = C(t)*Qgills*AFgills + Cfood*IRfood*AFfood - kelim*Cfish*Mfish
    """
    # Step 3: Calculate the three components of the accumulation equation
    uptake_gills = C_t * Qgills * AFgills
    uptake_food = Cfood * IRfood * AFfood
    elimination = kelim * Cfish * Mfish
    
    accumulation_rate = uptake_gills + uptake_food - elimination
    return accumulation_rate

# --- Calculations for PFOS ---
C_t_pfos = calculate_water_concentration(Cin_pfos, half_life_pfos_years, log_Kow_pfos)
accumulation_pfos = calculate_accumulation_rate(C_t_pfos, kelim_pfos)

# --- Calculations for PFOA ---
C_t_pfoa = calculate_water_concentration(Cin_pfoa, half_life_pfoa_years, log_Kow_pfoa)
accumulation_pfoa = calculate_accumulation_rate(C_t_pfoa, kelim_pfoa)

# --- Total Accumulation ---
total_accumulation = accumulation_pfos + accumulation_pfoa

# --- Final Output ---
print("This script calculates the net chemical accumulation rate in fish after 365 days.")
print("\n--- Individual Chemical Results ---")
print(f"PFOS accumulation rate: {accumulation_pfos:.2f} ng/day")
print(f"PFOA accumulation rate: {accumulation_pfoa:.2f} ng/day")
print("\n--- Final Equation ---")
print(f"Total Accumulation Rate = PFOS Rate + PFOA Rate")
print(f"Total Accumulation Rate = {accumulation_pfos:.2f} ng/day + {accumulation_pfoa:.2f} ng/day")
print(f"Total Accumulation Rate = {total_accumulation:.2f} ng/day")

# The final answer is the total accumulation rate.
# <<<6078314.93>>>
print("\nFinal Answer:")
print("The total chemical accumulation rate in the fish after 365 days is:")
print(f"{total_accumulation:.2f} ng/day")
print("<<<6078314.93>>>")