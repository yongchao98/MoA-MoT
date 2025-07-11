import math

def solve_environmental_problem():
    """
    Calculates the concentration of PFOS and PFOA in water and their 
    net accumulation rate in fish after 365 days.
    """
    # --- Part 1: Environmental and Fish Parameters ---

    # Environment
    V = 10000  # Volume of water body (L)
    Qin = 900  # Inflow rate (L/d)
    Qout = 1600  # Outflow rate (L/d)
    foc = 0.001  # Fraction of organic carbon (0.1%)
    t = 365  # time (days)

    # Fish
    Mfish = 1000  # Mass of fish (g)
    Cfood = 100  # Concentration in food (ng/g)
    IRfood = 20  # Ingestion rate of food (g/day)
    AFgills = 0.8  # Absorption fraction through gills
    AFfood = 0.9  # Absorption fraction for food
    Qgills = 100  # Gill flow rate (L/day)
    Cfish = 10  # Concentration in fish (ng/g)

    # Chemical-specific parameters
    params = {
        'PFOS': {
            'Cin': 2.6,        # ng/L
            'half_life': 91,   # years
            'log_Kow': 4.0,
            'kelim': 0.069     # days^-1
        },
        'PFOA': {
            'Cin': 211300,     # ng/L
            'half_life': 238,  # years
            'log_Kow': 4.5,
            'kelim': 0.023     # days^-1
        }
    }

    # --- Part 2: Calculations ---
    
    results = {}

    for chemical, p in params.items():
        # --- Calculate water concentration C(t) ---
        
        # Convert half-life to days
        half_life_days = p['half_life'] * 365.25
        
        # Degradation rate constant k (days^-1)
        k = math.log(2) / half_life_days
        
        # Hydraulic residence rate r (days^-1)
        r = Qout / V
        
        # Organic carbon partition coefficient Koc
        log_Koc = 0.81 * p['log_Kow'] + 0.01
        Koc = 10**log_Koc
        Kd = Koc # Assumption based on the problem context
        
        # Concentration in water C(t) at t=365 days, assuming C(0)=0
        # C(t) = (Cin * Qin / (Qout + Qin * (1 + Kd * foc))) * (1 - e^-(k + r)t)
        numerator = p['Cin'] * Qin
        denominator = Qout + Qin * (1 + Kd * foc)
        rate_sum = k + r
        exp_term = 1 - math.exp(-rate_sum * t)
        
        C_water_t = (numerator / denominator) * exp_term

        # --- Calculate POP accumulation rate in fish (ng/day) ---
        # Accumulation = Gill_Uptake + Food_Uptake - Elimination
        
        uptake_gills = C_water_t * Qgills * AFgills
        uptake_food = Cfood * IRfood * AFfood
        elimination = p['kelim'] * Cfish * Mfish
        
        net_accumulation_rate = uptake_gills + uptake_food - elimination
        
        results[chemical] = {
            "uptake_gills": uptake_gills,
            "uptake_food": uptake_food,
            "elimination": elimination,
            "net_accumulation_rate": net_accumulation_rate
        }

    # --- Part 3: Output the Results ---

    print("--- Chemical Accumulation Rate in Fish after 365 Days ---\n")

    # PFOS Results
    pfos_res = results['PFOS']
    print("For PFOS:")
    print(f"Net Accumulation Rate = (Uptake from Water) + (Uptake from Food) - (Elimination)")
    print(f"Net Accumulation Rate = {pfos_res['uptake_gills']:.2f} ng/day + {pfos_res['uptake_food']:.2f} ng/day - {pfos_res['elimination']:.2f} ng/day")
    print(f"Total Net Accumulation Rate for PFOS = {pfos_res['net_accumulation_rate']:.2f} ng/day\n")

    # PFOA Results
    pfoa_res = results['PFOA']
    print("For PFOA:")
    print(f"Net Accumulation Rate = (Uptake from Water) + (Uptake from Food) - (Elimination)")
    print(f"Net Accumulation Rate = {pfoa_res['uptake_gills']:.2f} ng/day + {pfoa_res['uptake_food']:.2f} ng/day - {pfoa_res['elimination']:.2f} ng/day")
    print(f"Total Net Accumulation Rate for PFOA = {pfoa_res['net_accumulation_rate']:.2f} ng/day")

solve_environmental_problem()