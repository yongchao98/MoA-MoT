import math

def solve_fish_accumulation():
    """
    Calculates the total chemical accumulation rate in fish for PFOS and PFOA.
    """
    # --- General Parameters ---
    # Environment
    V_env = 10000  # L
    Qin = 900  # L/d
    Qout = 1600  # L/d
    foc = 0.001  # 0.1% organic carbon fraction

    # Fish
    M_fish = 1000  # g
    C_food = 100  # ng/g
    IR_food = 20  # g/day
    AF_gills = 0.8
    AF_food = 0.9
    Q_gills = 100  # L/day

    # Time
    t = 365  # days
    C0_water = 0 # ng/L, assuming pristine environment initially

    # --- PFOS Calculation ---
    print("--- Calculating for PFOS ---")

    # PFOS specific parameters
    Cin_pfos = 2.6  # ng/L
    t_half_pfos_yr = 91  # years
    log_Kow_pfos = 4.0
    kelim_pfos = 0.069  # days^-1
    C_fish_pfos = 10  # ng/g

    # Step 1: Calculate water concentration C(t) for PFOS
    log_Koc_pfos = 0.81 * log_Kow_pfos + 0.01
    Koc_pfos = 10**log_Koc_pfos
    Kd_pfos = Koc_pfos * foc
    
    t_half_pfos_days = t_half_pfos_yr * 365
    k_pfos = math.log(2) / t_half_pfos_days
    r_env = Qout / V_env
    
    # Using the provided equation for C(t)
    # C(t) = (Cin * Qin / (Qout + Qin * (1 + Kd * foc))) * (1 - e^-(k + r)t)
    C_t_pfos_numerator = Cin_pfos * Qin
    C_t_pfos_denominator = Qout + Qin * (1 + Kd_pfos * foc)
    C_t_pfos_loss_rate = k_pfos + r_env
    C_t_pfos = (C_t_pfos_numerator / C_t_pfos_denominator) * (1 - math.exp(-C_t_pfos_loss_rate * t)) + C0_water * math.exp(-C_t_pfos_loss_rate * t)

    print("1. Water Concentration Equation for PFOS at t=365 days:")
    print(f"C_pfos(t) = ({Cin_pfos} * {Qin} / ({Qout} + {Qin} * (1 + {Kd_pfos:.3f} * {foc}))) * (1 - e^-(({k_pfos:.6f} + {r_env:.3f}) * {t}))")
    print(f"Result: C_pfos(365) = {C_t_pfos:.4f} ng/L\n")

    # Step 2: Calculate accumulation rate for PFOS
    # POP accumulation = C(t) * Qgills * AFgills + Cfood * IRfood * AFfood - kelim * Cfish * Mfish
    uptake_gills_pfos = C_t_pfos * Q_gills * AF_gills
    uptake_food_pfos = C_food * IR_food * AF_food
    elimination_pfos = kelim_pfos * C_fish_pfos * M_fish
    accum_rate_pfos = uptake_gills_pfos + uptake_food_pfos - elimination_pfos
    
    print("2. Accumulation Rate Equation for PFOS at t=365 days:")
    print(f"Accum_pfos = {C_t_pfos:.4f} * {Q_gills} * {AF_gills} + {C_food} * {IR_food} * {AF_food} - {kelim_pfos} * {C_fish_pfos} * {M_fish}")
    print(f"Result: Accum_pfos = {accum_rate_pfos:.4f} ng/day\n")

    # --- PFOA Calculation ---
    print("--- Calculating for PFOA ---")

    # PFOA specific parameters
    Cin_pfoa = 211300  # ng/L
    t_half_pfoa_yr = 238  # years
    log_Kow_pfoa = 4.5
    kelim_pfoa = 0.023  # days^-1
    C_fish_pfoa = 10  # ng/g

    # Step 1: Calculate water concentration C(t) for PFOA
    log_Koc_pfoa = 0.81 * log_Kow_pfoa + 0.01
    Koc_pfoa = 10**log_Koc_pfoa
    Kd_pfoa = Koc_pfoa * foc
    
    t_half_pfoa_days = t_half_pfoa_yr * 365
    k_pfoa = math.log(2) / t_half_pfoa_days
    # r_env is the same for both chemicals
    
    # Using the provided equation for C(t)
    C_t_pfoa_numerator = Cin_pfoa * Qin
    C_t_pfoa_denominator = Qout + Qin * (1 + Kd_pfoa * foc)
    C_t_pfoa_loss_rate = k_pfoa + r_env
    C_t_pfoa = (C_t_pfoa_numerator / C_t_pfoa_denominator) * (1 - math.exp(-C_t_pfoa_loss_rate * t)) + C0_water * math.exp(-C_t_pfoa_loss_rate * t)

    print("1. Water Concentration Equation for PFOA at t=365 days:")
    print(f"C_pfoa(t) = ({Cin_pfoa} * {Qin} / ({Qout} + {Qin} * (1 + {Kd_pfoa:.3f} * {foc}))) * (1 - e^-(({k_pfoa:.6f} + {r_env:.3f}) * {t}))")
    print(f"Result: C_pfoa(365) = {C_t_pfoa:.4f} ng/L\n")

    # Step 2: Calculate accumulation rate for PFOA
    uptake_gills_pfoa = C_t_pfoa * Q_gills * AF_gills
    uptake_food_pfoa = C_food * IR_food * AF_food
    elimination_pfoa = kelim_pfoa * C_fish_pfoa * M_fish
    accum_rate_pfoa = uptake_gills_pfoa + uptake_food_pfoa - elimination_pfoa
    
    print("2. Accumulation Rate Equation for PFOA at t=365 days:")
    print(f"Accum_pfoa = {C_t_pfoa:.4f} * {Q_gills} * {AF_gills} + {C_food} * {IR_food} * {AF_food} - {kelim_pfoa} * {C_fish_pfoa} * {M_fish}")
    print(f"Result: Accum_pfoa = {accum_rate_pfoa:.4f} ng/day\n")
    
    # --- Total Calculation ---
    total_accum_rate = accum_rate_pfos + accum_rate_pfoa
    
    print("--- Total Accumulation Rate ---")
    print(f"Total Net Accumulation Rate = Accum_pfos + Accum_pfoa")
    print(f"Total Net Accumulation Rate = {accum_rate_pfos:.4f} ng/day + {accum_rate_pfoa:.4f} ng/day")
    print(f"Total Net Accumulation Rate = {total_accum_rate:.4f} ng/day")
    
    return total_accum_rate

# Run the calculation and print the final answer in the required format
final_answer = solve_fish_accumulation()
print(f"\nThe final combined chemical accumulation rate in the fish at 365 days is {final_answer:.2f} ng/day.")
print(f"<<<{final_answer:.2f}>>>")
