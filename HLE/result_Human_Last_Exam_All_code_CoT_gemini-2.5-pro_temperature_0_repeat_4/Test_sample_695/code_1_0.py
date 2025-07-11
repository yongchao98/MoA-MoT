import math

def solve_accumulation():
    """
    Calculates the net accumulation rate of PFOS and PFOA in fish after 365 days.
    """
    # --- Given Constants ---

    # Environment
    V_water = 10000  # L
    Qin = 900      # L/d
    Qout = 1600    # L/d
    t = 365        # days

    # Fish
    M_fish = 1000    # g
    C_food = 100     # ng/g (for each chemical)
    IR_food = 20     # g/day
    AF_gills = 0.8
    AF_food = 0.9
    Q_gills = 100    # L/day
    C_fish_initial = 10 # ng/g

    # PFOS Data
    Cin_pfos = 2.6      # ng/L
    thalf_pfos_yr = 91  # years
    kelim_pfos = 0.069  # days^-1

    # PFOA Data
    Cin_pfoa = 211300   # ng/L
    thalf_pfoa_yr = 238 # years
    kelim_pfoa = 0.023  # days^-1

    # --- Step 1: Calculate Degradation Rate Constants (k) ---
    thalf_pfos_days = thalf_pfos_yr * 365.25
    k_pfos = math.log(2) / thalf_pfos_days

    thalf_pfoa_days = thalf_pfoa_yr * 365.25
    k_pfoa = math.log(2) / thalf_pfoa_days

    # --- Step 2: Calculate Steady-State Water Concentration (C_water_ss) ---
    # Using standard CSTR model: Css = (Cin * Qin) / (Qout + k * V)
    C_water_pfos = (Cin_pfos * Qin) / (Qout + k_pfos * V_water)
    C_water_pfoa = (Cin_pfoa * Qin) / (Qout + k_pfoa * V_water)

    # --- Step 3: Calculate Fish Concentration at t=365 days (C_fish) ---
    # Solves the differential equation for bioaccumulation over time.
    
    # For PFOS
    uptake_rate_pfos = (C_water_pfos * Q_gills * AF_gills) + (C_food * IR_food * AF_food)
    C_fish_ss_pfos = uptake_rate_pfos / (M_fish * kelim_pfos)
    C_fish_365_pfos = C_fish_ss_pfos * (1 - math.exp(-kelim_pfos * t)) + C_fish_initial * math.exp(-kelim_pfos * t)

    # For PFOA
    uptake_rate_pfoa = (C_water_pfoa * Q_gills * AF_gills) + (C_food * IR_food * AF_food)
    C_fish_ss_pfoa = uptake_rate_pfoa / (M_fish * kelim_pfoa)
    C_fish_365_pfoa = C_fish_ss_pfoa * (1 - math.exp(-kelim_pfoa * t)) + C_fish_initial * math.exp(-kelim_pfoa * t)

    # --- Step 4: Calculate Final Net Accumulation Rate at t=365 days ---
    # Accumulation = Uptake - Elimination
    
    # PFOS Calculation
    uptake_gills_pfos = C_water_pfos * Q_gills * AF_gills
    uptake_food_pfos = C_food * IR_food * AF_food
    elimination_pfos = kelim_pfos * C_fish_365_pfos * M_fish
    accumulation_pfos = uptake_gills_pfos + uptake_food_pfos - elimination_pfos

    # PFOA Calculation
    uptake_gills_pfoa = C_water_pfoa * Q_gills * AF_gills
    uptake_food_pfoa = C_food * IR_food * AF_food
    elimination_pfoa = kelim_pfoa * C_fish_365_pfoa * M_fish
    accumulation_pfoa = uptake_gills_pfoa + uptake_food_pfoa - elimination_pfoa

    # --- Print Results ---
    print("Based on the calculations, the net accumulation rates at day 365 are:\n")
    
    print("PFOS Accumulation Rate Calculation:")
    print(f"= (C_water * Q_gills * AF_gills) + (C_food * IR_food * AF_food) - (kelim * C_fish_365 * M_fish)")
    print(f"= ({C_water_pfos:.3f} ng/L * {Q_gills} L/day * {AF_gills}) + ({C_food} ng/g * {IR_food} g/day * {AF_food}) - ({kelim_pfos} 1/day * {C_fish_365_pfos:.3f} ng/g * {M_fish} g)")
    print(f"= ({uptake_gills_pfos:.3f} ng/day) + ({uptake_food_pfos:.3f} ng/day) - ({elimination_pfos:.3f} ng/day)")
    print(f"= {accumulation_pfos:.3f} ng/day\n")

    print("PFOA Accumulation Rate Calculation:")
    print(f"= (C_water * Q_gills * AF_gills) + (C_food * IR_food * AF_food) - (kelim * C_fish_365 * M_fish)")
    print(f"= ({C_water_pfoa:.3f} ng/L * {Q_gills} L/day * {AF_gills}) + ({C_food} ng/g * {IR_food} g/day * {AF_food}) - ({kelim_pfoa} 1/day * {C_fish_365_pfoa:.3f} ng/g * {M_fish} g)")
    print(f"= ({uptake_gills_pfoa:.3f} ng/day) + ({uptake_food_pfoa:.3f} ng/day) - ({elimination_pfoa:.3f} ng/day)")
    print(f"= {accumulation_pfoa:.3f} ng/day")


solve_accumulation()