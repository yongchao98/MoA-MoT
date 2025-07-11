import math

# --- Given Parameters ---

# Environment
V = 10000  # L (Volume of freshwater body)
Qin = 900  # L/d (Inflow rate)
Qout = 1600 # L/d (Outflow rate)
foc = 0.001 # 0.1% organic carbon
C0_water = 0 # ng/L (Initial concentration in water, assumed pristine)
t = 365 # days

# Fish
Mfish = 1000 # g (Fish weight)
Cfood = 100 # ng/g (Chemical concentration in food)
IRfood = 20 # g/day (Ingestion rate)
AFgills = 0.8 # Absorption fraction gills
AFfood = 0.9 # Absorption fraction food
Qgills = 100 # L/day (Gill flow rate)
Cfish_initial = 10 # ng/g (Initial concentration in fish)

# --- PFOS Parameters ---
Cin_pfos = 2.6 # ng/L
half_life_pfos_years = 91
logKow_pfos = 4.0
kelim_pfos = 0.069 # days^-1

# --- PFOA Parameters ---
Cin_pfoa = 211300 # ng/L
half_life_pfoa_years = 238
logKow_pfoa = 4.5
kelim_pfoa = 0.023 # days^-1


print("Step 1: Calculate Water Concentration C(t=365) for PFOS and PFOA\n")

# --- Calculations for PFOS ---
print("--- For PFOS ---")
# Water concentration
half_life_pfos_days = half_life_pfos_years * 365
k_decay_pfos = math.log(2) / half_life_pfos_days
logKoc_pfos = 0.81 * logKow_pfos + 0.01
Koc_pfos = 10**logKoc_pfos
Kd_pfos = Koc_pfos * foc
r_hydraulic = Qout / V
k_plus_r_pfos = k_decay_pfos + r_hydraulic
denom_pfos = Qout + Qin * (1 + Kd_pfos * foc)
exp_term_pfos = math.exp(-k_plus_r_pfos * t)
C_water_365_pfos = (Cin_pfos * Qin / denom_pfos) * (1 - exp_term_pfos) + C0_water * exp_term_pfos
print(f"Water concentration of PFOS, C_water(365):")
print(f"Equation: C(t) = (Cin * Qin / (Qout + Qin * (1 + Kd * foc))) * (1 - e^-(k + r)t) + C0 * e^-(k + r)t")
print(f"C(365) = ({Cin_pfos:.1f} * {Qin} / ({Qout} + {Qin} * (1 + {Kd_pfos:.3f} * {foc:.3f}))) * (1 - e^-({k_plus_r_pfos:.6f} * {t})) + {C0_water} * ...")
print(f"C(365) = {C_water_365_pfos:.4f} ng/L\n")

# Fish accumulation
M_pop_initial_pfos = Cfish_initial * Mfish
uptake_gills_pfos = C_water_365_pfos * Qgills * AFgills
uptake_food_pfos = Cfood * IRfood * AFfood
U_pfos = uptake_gills_pfos + uptake_food_pfos
M_pop_ss_pfos = U_pfos / kelim_pfos # Steady-state mass
exp_term_fish_pfos = math.exp(-kelim_pfos * t)
M_pop_365_pfos = M_pop_ss_pfos * (1 - exp_term_fish_pfos) + M_pop_initial_pfos * exp_term_fish_pfos
C_fish_365_pfos = M_pop_365_pfos / Mfish

print(f"Final concentration of PFOS in fish, C_fish(365):")
print(f"Equation: Mass_pop(t) = (U / kelim) * (1 - e^(-kelim * t)) + Mass_pop(0) * e^(-kelim * t)")
print(f"U_pfos = {C_water_365_pfos:.4f} * {Qgills} * {AFgills} + {Cfood} * {IRfood} * {AFfood} = {U_pfos:.4f} ng/day")
print(f"Mass_pop(365) = ({U_pfos:.4f} / {kelim_pfos}) * (1 - e^(-{kelim_pfos} * {t})) + {M_pop_initial_pfos} * e^(-{kelim_pfos} * {t}) = {M_pop_365_pfos:.4f} ng")
print(f"C_fish(365) = {M_pop_365_pfos:.4f} ng / {Mfish} g = {C_fish_365_pfos:.4f} ng/g\n")

# --- Calculations for PFOA ---
print("--- For PFOA ---")
# Water concentration
half_life_pfoa_days = half_life_pfoa_years * 365
k_decay_pfoa = math.log(2) / half_life_pfoa_days
logKoc_pfoa = 0.81 * logKow_pfoa + 0.01
Koc_pfoa = 10**logKoc_pfoa
Kd_pfoa = Koc_pfoa * foc
k_plus_r_pfoa = k_decay_pfoa + r_hydraulic
denom_pfoa = Qout + Qin * (1 + Kd_pfoa * foc)
exp_term_pfoa = math.exp(-k_plus_r_pfoa * t)
C_water_365_pfoa = (Cin_pfoa * Qin / denom_pfoa) * (1 - exp_term_pfoa) + C0_water * exp_term_pfoa
print(f"Water concentration of PFOA, C_water(365):")
print(f"Equation: C(t) = (Cin * Qin / (Qout + Qin * (1 + Kd * foc))) * (1 - e^-(k + r)t) + C0 * e^-(k + r)t")
print(f"C(365) = ({Cin_pfoa} * {Qin} / ({Qout} + {Qin} * (1 + {Kd_pfoa:.3f} * {foc:.3f}))) * (1 - e^-({k_plus_r_pfoa:.6f} * {t})) + {C0_water} * ...")
print(f"C(365) = {C_water_365_pfoa:.4f} ng/L\n")

# Fish accumulation
M_pop_initial_pfoa = Cfish_initial * Mfish
uptake_gills_pfoa = C_water_365_pfoa * Qgills * AFgills
uptake_food_pfoa = Cfood * IRfood * AFfood
U_pfoa = uptake_gills_pfoa + uptake_food_pfoa
M_pop_ss_pfoa = U_pfoa / kelim_pfoa # Steady-state mass
exp_term_fish_pfoa = math.exp(-kelim_pfoa * t)
M_pop_365_pfoa = M_pop_ss_pfoa * (1 - exp_term_fish_pfoa) + M_pop_initial_pfoa * exp_term_fish_pfoa
C_fish_365_pfoa = M_pop_365_pfoa / Mfish

print(f"Final concentration of PFOA in fish, C_fish(365):")
print(f"Equation: Mass_pop(t) = (U / kelim) * (1 - e^(-kelim * t)) + Mass_pop(0) * e^(-kelim * t)")
print(f"U_pfoa = {C_water_365_pfoa:.4f} * {Qgills} * {AFgills} + {Cfood} * {IRfood} * {AFfood} = {U_pfoa:.4f} ng/day")
print(f"Mass_pop(365) = ({U_pfoa:.4f} / {kelim_pfoa}) * (1 - e^(-{kelim_pfoa} * {t})) + {M_pop_initial_pfoa} * e^(-{kelim_pfoa} * {t}) = {M_pop_365_pfoa:.4f} ng")
print(f"C_fish(365) = {M_pop_365_pfoa:.4f} ng / {Mfish} g = {C_fish_365_pfoa:.4f} ng/g\n")


# --- Final Summation ---
total_concentration = C_fish_365_pfos + C_fish_365_pfoa
print("--- Total Accumulated Concentration ---")
print(f"The total accumulated concentration is the sum of the final concentrations of PFOS and PFOA.")
print(f"Total Concentration = C_fish_PFOS(365) + C_fish_PFOA(365)")
print(f"Total Concentration = {C_fish_365_pfos:.4f} ng/g + {C_fish_365_pfoa:.4f} ng/g = {total_concentration:.4f} ng/g")
print(f"\nFinal Answer: {total_concentration:.1f}")