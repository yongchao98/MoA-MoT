import math

def calculate_fish_accumulation():
    """
    Calculates the accumulation of PFOS and PFOA in fish based on provided environmental models.
    """
    # Part 1: Define all given parameters
    # Environment parameters
    V = 10000.0      # L (Volume)
    Qin = 900.0      # L/d (Inflow rate)
    Qout = 1600.0    # L/d (Outflow rate)
    foc = 0.001      # unitless (0.1% organic carbon)
    t = 365.0        # days
    C0 = 0.0         # ng/L (Initial concentration in pristine environment)

    # Fish parameters
    Mfish = 1000.0   # g (Weight)
    Cfood = 100.0    # ng/g (Food concentration)
    IRfood = 20.0    # g/day (Ingestion rate)
    AFgills = 0.8    # unitless (Absorption fraction, gills)
    AFfood = 0.9     # unitless (Absorption fraction, food)
    Qgills = 100.0   # L/day (Gill flow rate)
    Cfish = 10.0     # ng/g (Current fish concentration)

    # Chemical-specific parameters
    chemicals = {
        'PFOS': {
            'Cin': 2.6,        # ng/L
            't_half_years': 91,
            'log_Kow': 4.0,
            'kelim': 0.069   # days⁻¹
        },
        'PFOA': {
            'Cin': 211300.0,   # ng/L
            't_half_years': 238,
            'log_Kow': 4.5,
            'kelim': 0.023   # days⁻¹
        }
    }

    results = {}
    print("This script calculates the net accumulation rate (in ng/day) of PFOS and PFOA in fish after 365 days.\n")

    # Part 2: Loop through each chemical to perform calculations
    for name, params in chemicals.items():
        print(f"=============== CALCULATING FOR {name} ===============\n")

        # Step A: Calculate water concentration C(t) at t=365 days.
        print("Step 1: Calculate Water Concentration C(t)")
        
        # A.1: Calculate log Koc and Koc
        log_Koc = 0.81 * params['log_Kow'] + 0.01
        Koc = 10**log_Koc
        print(f"log Koc = 0.81 * {params['log_Kow']} + 0.01 = {log_Koc:.4f}")
        print(f"Koc = 10^{log_Koc:.4f} = {Koc:.4f}")

        # A.2: Calculate Kd
        Kd = Koc * foc
        print(f"Kd = Koc * foc = {Koc:.4f} * {foc} = {Kd:.4f}")

        # A.3: Calculate degradation rate constant k
        t_half_days = params['t_half_years'] * 365.0
        k = math.log(2) / t_half_days
        print(f"k = ln(2) / ({params['t_half_years']} * 365 days) = {k:.4e} days⁻¹")

        # A.4: Calculate r using the provided formula
        denominator_term = Qout + Qin * (1 + Kd * foc)
        r = denominator_term / V
        print(f"r = (Qout + Qin * (1 + Kd * foc)) / V")
        print(f"r = ({Qout} + {Qin} * (1 + {Kd:.4f} * {foc})) / {V} = {r:.4f} days⁻¹")

        # A.5: Calculate C(t)
        # The term (Cin * Qin / denominator_term) is the steady-state concentration
        Css = (params['Cin'] * Qin) / denominator_term
        exponent = - (k + r) * t
        Ct = Css * (1 - math.exp(exponent)) + C0 * math.exp(exponent)
        
        print("\nUsing equation: C(t) = (Cin*Qin / (Qout+Qin*(1+Kd*foc))) * (1-e^-(k+r)t) + C0*e^-(k+r)t")
        print(f"C(365) = ({params['Cin']} * {Qin} / ({Qout} + {Qin} * (1 + {Kd:.4f} * {foc}))) * (1-e^-(({k:.4e})+({r:.4f}))*{t})")
        print(f"Water Concentration C(365) for {name} = {Ct:.4f} ng/L\n")

        # Step B: Calculate POP accumulation rate in fish
        print("Step 2: Calculate POP Accumulation Rate in Fish")
        uptake_gills = Ct * Qgills * AFgills
        uptake_food = Cfood * IRfood * AFfood
        elimination = params['kelim'] * Cfish * Mfish
        
        pop_accumulation = uptake_gills + uptake_food - elimination
        results[name] = pop_accumulation
        
        print("Using equation: POP accumulation = C(t)*Qgills*AFgills + Cfood*IRfood*AFfood - kelim*Cfish*Mfish")
        print(f"POP accumulation = {Ct:.4f}*{Qgills}*{AFgills} + {Cfood}*{IRfood}*{AFfood} - {params['kelim']}*{Cfish}*{Mfish}")
        print(f"POP accumulation = {uptake_gills:.2f} (from gills) + {uptake_food:.2f} (from food) - {elimination:.2f} (elimination)")
        print(f"\nFinal Net Accumulation Rate for {name}: {pop_accumulation:.2f} ng/day\n")

    return results

# Run the calculation
final_results = calculate_fish_accumulation()
pfos_accumulation_rate = final_results['PFOS']
pfoa_accumulation_rate = final_results['PFOA']
<<<PFOS Rate: 1184.83 ng/day, PFOA Rate: 6077130.29 ng/day>>>