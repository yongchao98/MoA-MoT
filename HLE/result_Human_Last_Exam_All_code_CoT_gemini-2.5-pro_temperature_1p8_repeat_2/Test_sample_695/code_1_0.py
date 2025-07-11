import math

def solve_environmental_problem():
    # --- Input Parameters ---

    # Environment
    V_water = 10000  # L
    Qin = 900  # L/d
    Qout = 1600  # L/d
    foc = 0.001  # 0.1% organic carbon
    C0_water = 0 # ng/L (pristine environment)
    t_days = 365 # days

    # PFOS Parameters
    Cin_pfos = 2.6  # ng/L
    t_half_pfos = 91 * 365  # days
    logKow_pfos = 4.0
    kelim_pfos = 0.069  # days⁻¹

    # PFOA Parameters
    Cin_pfoa = 211300  # ng/L
    t_half_pfoa = 238 * 365  # days
    logKow_pfoa = 4.5
    kelim_pfoa = 0.023  # days⁻¹

    # Fish Parameters
    M_fish = 1000  # g
    C_food = 100  # ng/g (for each chemical)
    IR_food = 20  # g/day
    AF_gills = 0.8
    AF_food = 0.9
    Q_gills = 100  # L/day
    C_fish_initial = 10  # ng/g

    # --- Part 1: Calculate Water Concentration C(t) at t=365 days ---

    # Shared calculation for water
    r = Qout / V_water  # d⁻¹

    # Calculations for PFOS in water
    k_pfos = math.log(2) / t_half_pfos
    logKoc_pfos = 0.81 * logKow_pfos + 0.01
    Koc_pfos = 10**logKoc_pfos
    Kd_pfos = Koc_pfos # Assume Kd = Koc
    denom_pfos = Qout + Qin * (1 + Kd_pfos * foc)
    C_water_ss_pfos = (Cin_pfos * Qin) / denom_pfos
    exponent_pfos = -(k_pfos + r) * t_days
    C_water_365_pfos = C_water_ss_pfos * (1 - math.exp(exponent_pfos)) + C0_water * math.exp(exponent_pfos)

    # Calculations for PFOA in water
    k_pfoa = math.log(2) / t_half_pfoa
    logKoc_pfoa = 0.81 * logKow_pfoa + 0.01
    Koc_pfoa = 10**logKoc_pfoa
    Kd_pfoa = Koc_pfoa # Assume Kd = Koc
    denom_pfoa = Qout + Qin * (1 + Kd_pfoa * foc)
    C_water_ss_pfoa = (Cin_pfoa * Qin) / denom_pfoa
    exponent_pfoa = -(k_pfoa + r) * t_days
    C_water_365_pfoa = C_water_ss_pfoa * (1 - math.exp(exponent_pfoa)) + C0_water * math.exp(exponent_pfoa)

    # --- Part 2: Calculate Fish Concentration C_fish(t) at t=365 days ---

    Mass_fish_initial = C_fish_initial * M_fish

    # Calculations for PFOS in fish
    Uptake_pfos = C_water_365_pfos * Q_gills * AF_gills + C_food * IR_food * AF_food
    Mass_ss_pfos = Uptake_pfos / kelim_pfos
    exponent_fish_pfos = -kelim_pfos * t_days
    Mass_365_pfos = Mass_ss_pfos * (1 - math.exp(exponent_fish_pfos)) + Mass_fish_initial * math.exp(exponent_fish_pfos)
    C_fish_365_pfos = Mass_365_pfos / M_fish

    # Calculations for PFOA in fish
    Uptake_pfoa = C_water_365_pfoa * Q_gills * AF_gills + C_food * IR_food * AF_food
    Mass_ss_pfoa = Uptake_pfoa / kelim_pfoa
    exponent_fish_pfoa = -kelim_pfoa * t_days
    Mass_365_pfoa = Mass_ss_pfoa * (1 - math.exp(exponent_fish_pfoa)) + Mass_fish_initial * math.exp(exponent_fish_pfoa)
    C_fish_365_pfoa = Mass_365_pfoa / M_fish
    
    # --- Final Output ---
    
    print("--- Analysis for PFOS at 365 days ---")
    print("\nFinal Concentration in Fish:")
    print(f"The accumulated concentration of PFOS in fish is {C_fish_365_pfos:.2f} ng/g.")
    
    print("\nAccumulation Rate Equation at day 365:")
    print("POP accumulation = C(t) * Qgills * AFgills + Cfood * IRfood * AFfood - kelim * Cfish * Mfish")
    elimination_term_pfos = kelim_pfos * C_fish_365_pfos * M_fish
    accumulation_rate_pfos = Uptake_pfos - elimination_term_pfos
    print(f"Rate (ng/day) = {C_water_365_pfos:.3f} ng/L * {Q_gills} L/day * {AF_gills} + {C_food} ng/g * {IR_food} g/day * {AF_food} - {kelim_pfos} day⁻¹ * {C_fish_365_pfos:.2f} ng/g * {M_fish} g")
    print(f"Rate (ng/day) = {accumulation_rate_pfos:.2f} ng/day")

    print("\n" + "="*40 + "\n")

    print("--- Analysis for PFOA at 365 days ---")
    print("\nFinal Concentration in Fish:")
    print(f"The accumulated concentration of PFOA in fish is {C_fish_365_pfoa:.2f} ng/g.")
    
    print("\nAccumulation Rate Equation at day 365:")
    print("POP accumulation = C(t) * Qgills * AFgills + Cfood * IRfood * AFfood - kelim * Cfish * Mfish")
    elimination_term_pfoa = kelim_pfoa * C_fish_365_pfoa * M_fish
    accumulation_rate_pfoa = Uptake_pfoa - elimination_term_pfoa
    print(f"Rate (ng/day) = {C_water_365_pfoa:.3f} ng/L * {Q_gills} L/day * {AF_gills} + {C_food} ng/g * {IR_food} g/day * {AF_food} - {kelim_pfoa} day⁻¹ * {C_fish_365_pfoa:.2f} ng/g * {M_fish} g")
    print(f"Rate (ng/day) = {accumulation_rate_pfoa:.2f} ng/day")

solve_environmental_problem()
<<<The final concentration of PFOS is 26.75 ng/g and the final concentration of PFOA is 100927.46 ng/g.>>>