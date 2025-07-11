import math

# Step 1: Define all given parameters

# Environment
V = 10000  # L (Volume)
Qin = 900  # L/d (Inflow rate)
Qout = 1600 # L/d (Outflow rate)
foc = 0.001 # 0.1% organic carbon
t = 365  # days
C0 = 0   # ng/L (Initial concentration in pristine environment)

# Fish parameters (common for both chemicals)
Mfish = 1000   # g (Fish weight)
Cfood = 100    # ng/g (Concentration in food)
IRfood = 20    # g/day (Ingestion rate)
AFgills = 0.8  # (Absorption fraction gills)
AFfood = 0.9   # (Absorption fraction food)
Qgills = 100   # L/day (Gill flow rate)
Cfish = 10     # ng/g (Initial concentration in fish)

# PFOS parameters
Cin_pfos = 2.6       # ng/L
half_life_pfos = 91  # years
log_kow_pfos = 4.0
kelim_pfos = 0.069   # days⁻¹

# PFOA parameters
Cin_pfoa = 211300      # ng/L
half_life_pfoa = 238 # years
log_kow_pfoa = 4.5
kelim_pfoa = 0.023   # days⁻¹


# Step 2: Function to calculate water concentration C(t)
def calculate_water_concentration(Cin, half_life_years, log_kow):
    # Convert half-life from years to days
    half_life_days = half_life_years * 365.25
    
    # Calculate degradation rate constant (k)
    k = math.log(2) / half_life_days
    
    # Calculate log Koc and then Koc
    log_koc = 0.81 * log_kow + 0.01
    Koc = 10**log_koc
    
    # Calculate partition coefficient (Kd)
    Kd = Koc * foc
    
    # Calculate residence rate (r)
    r = Qout / V
    
    # Calculate C(t) using the provided formula
    # C(t) = (Cin * Qin / (Qout + Qin * (1 + Kd * foc))) * (1 - e^-(k + r)t) + C0 * e^-(k + r)t
    # Since C0 = 0, the second term is zero.
    
    C_ss_term = (Cin * Qin) / (Qout + Qin * (1 + Kd * foc))
    exponent_term = math.exp(-(k + r) * t)
    Ct = C_ss_term * (1 - exponent_term) + C0 * exponent_term
    
    return Ct

# Step 3: Function to calculate POP accumulation rate in fish
def calculate_fish_accumulation_rate(Ct, kelim):
    uptake_gills = Ct * Qgills * AFgills
    uptake_food = Cfood * IRfood * AFfood
    elimination = kelim * Cfish * Mfish
    
    accumulation_rate = uptake_gills + uptake_food - elimination
    return accumulation_rate, uptake_gills, uptake_food, elimination

# Perform calculations for PFOS
Ct_pfos = calculate_water_concentration(Cin_pfos, half_life_pfos, log_kow_pfos)
acc_rate_pfos, ug_pfos, uf_pfos, e_pfos = calculate_fish_accumulation_rate(Ct_pfos, kelim_pfos)

# Perform calculations for PFOA
Ct_pfoa = calculate_water_concentration(Cin_pfoa, half_life_pfoa, log_kow_pfoa)
acc_rate_pfoa, ug_pfoa, uf_pfoa, e_pfoa = calculate_fish_accumulation_rate(Ct_pfoa, kelim_pfoa)

# Step 4: Sum the results and print the final equation
total_accumulation_rate = acc_rate_pfos + acc_rate_pfoa

print("This calculation determines the net rate of accumulation (in ng/day) for the pollutants in fish at the 365-day mark.\n")
print(f"Calculated water concentration for PFOS at 365 days (C(t)_pfos): {Ct_pfos:.4f} ng/L")
print(f"Calculated water concentration for PFOA at 365 days (C(t)_pfoa): {Ct_pfoa:.4f} ng/L\n")

print("The total accumulation rate is the sum of the rates for PFOS and PFOA.")
print("Total Accumulation Rate = (Accumulation Rate_PFOS) + (Accumulation Rate_PFOA)")
print("Accumulation Rate = C(t) * Qgills * AFgills + Cfood * IRfood * AFfood - kelim * Cfish * Mfish\n")

# Print the final detailed equation with all numbers
print("Final Equation:")
print(f"Total Accumulation (ng/day) = \n"
      f"  (C(t)_pfos * {Qgills} * {AFgills} + {Cfood} * {IRfood} * {AFfood} - {kelim_pfos} * {Cfish} * {Mfish}) + \n"
      f"  (C(t)_pfoa * {Qgills} * {AFgills} + {Cfood} * {IRfood} * {AFfood} - {kelim_pfoa} * {Cfish} * {Mfish})\n")

print("Plugging in the calculated values for C(t):")
print(f"Total Accumulation (ng/day) = \n"
      f"  ({Ct_pfos:.4f} * {Qgills} * {AFgills} + {Cfood} * {IRfood} * {AFfood} - {kelim_pfos} * {Cfish} * {Mfish}) + \n"
      f"  ({Ct_pfoa:.4f} * {Qgills} * {AFgills} + {Cfood} * {IRfood} * {AFfood} - {kelim_pfoa} * {Cfish} * {Mfish})\n")

print("Calculating each component:")
print(f"Total Accumulation (ng/day) = ({ug_pfos:.2f} + {uf_pfos:.2f} - {e_pfos:.2f}) + ({ug_pfoa:.2f} + {uf_pfoa:.2f} - {e_pfoa:.2f})")
print(f"Total Accumulation (ng/day) = ({acc_rate_pfos:.2f}) + ({acc_rate_pfoa:.2f})\n")

print(f"Final Result: {total_accumulation_rate:.2f} ng/day")
print(f"<<<{total_accumulation_rate:.2f}>>>")