import math

def calculate_chemical_accumulation():
    """
    This script calculates the net accumulation rate of PFOS and PFOA in fish
    based on the environmental and biological parameters provided.
    """

    # --- Step 1: Define all given parameters ---

    # Environment Parameters
    V = 10000       # Volume of freshwater body (L)
    Qin = 900       # Inflow rate (L/d)
    Qout = 1600     # Outflow rate (L/d)
    foc = 0.001     # Fraction of organic carbon (0.1%)
    t = 365         # Time (days)
    C0 = 0          # Initial water concentration (ng/L), assumed pristine

    # Fish Parameters
    Mfish = 1000    # Fish weight (g)
    Cfood = 100     # Food concentration (ng/g)
    IRfood = 20     # Ingestion rate (g/day)
    AFgills = 0.8   # Absorption fraction for gills
    AFfood = 0.9    # Absorption fraction for food
    Qgills = 100    # Gill flow rate (L/day)

    # Chemical-specific parameters in a list of dictionaries
    chemicals = [
        {
            "name": "PFOS",
            "Cin": 2.6,        # Inflow concentration (ng/L)
            "t_half_years": 91,
            "log_Kow": 4.0,
            "kelim": 0.069,    # Elimination rate constant (days⁻¹)
            "Cfish": 10        # Fish concentration (ng/g)
        },
        {
            "name": "PFOA",
            "Cin": 211300,
            "t_half_years": 238,
            "log_Kow": 4.5,
            "kelim": 0.023,
            "Cfish": 10
        }
    ]

    final_results = {}

    for params in chemicals:
        name = params["name"]
        Cin = params["Cin"]
        t_half_years = params["t_half_years"]
        log_Kow = params["log_Kow"]
        kelim = params["kelim"]
        Cfish = params["Cfish"]

        print(f"--- Calculations for {name} ---")

        # --- Step 2: Calculate intermediate values ---
        
        # Degradation rate constant (k)
        t_half_days = t_half_years * 365
        k = math.log(2) / t_half_days

        # Koc from log Koc
        log_Koc = 0.81 * log_Kow + 0.01
        Koc = 10**log_Koc

        # Distribution coefficient (Kd)
        Kd = Koc * foc

        # Hydraulic residence rate (r)
        r = Qout / V

        # --- Step 3: Calculate water concentration C(t) ---
        print("\nCalculating water concentration at 365 days, C(365):")
        print(f"Using equation: C(t) = (Cin * Qin / (Qout + Qin * (1 + Kd * foc))) * (1 - e^-(k + r)t)")
        
        k_plus_r = k + r
        exp_term = math.exp(-k_plus_r * t)
        denominator_Ct = Qout + Qin * (1 + Kd * foc)
        C_t = (Cin * Qin / denominator_Ct) * (1 - exp_term) + C0 * exp_term
        
        print(f"C(365) for {name} = ({Cin} * {Qin} / ({Qout} + {Qin} * (1 + {Kd:.4f} * {foc}))) * (1 - e^-({k:.6f} + {r})*{t})")
        print(f"Resulting water concentration C(365): {C_t:.4f} ng/L\n")

        # --- Step 4: Calculate net accumulation rate in fish ---
        print("Calculating net accumulation rate in fish at 365 days:")
        print("Using equation: Rate = C(t)*Qgills*AFgills + Cfood*IRfood*AFfood - kelim*Cfish*Mfish")

        uptake_gills = C_t * Qgills * AFgills
        uptake_food = Cfood * IRfood * AFfood
        elimination = kelim * Cfish * Mfish
        
        net_accumulation_rate = uptake_gills + uptake_food - elimination
        final_results[name] = net_accumulation_rate

        print(f"Uptake from gills = {C_t:.4f} * {Qgills} * {AFgills} = {uptake_gills:.2f} ng/day")
        print(f"Uptake from food = {Cfood} * {IRfood} * {AFfood} = {uptake_food:.2f} ng/day")
        print(f"Elimination = {kelim} * {Cfish} * {Mfish} = {elimination:.2f} ng/day")
        print(f"Final Net Accumulation Rate for {name} = {uptake_gills:.2f} + {uptake_food:.2f} - {elimination:.2f} = {net_accumulation_rate:.2f} ng/day")
        print("\n" + "="*50 + "\n")

    # --- Step 5: Present the final answers ---
    print("\n--- Summary of Final Answers ---")
    print(f"The net chemical accumulation rate for PFOS in fish after 365 days is {final_results['PFOS']:.2f} ng/day.")
    print(f"The net chemical accumulation rate for PFOA in fish after 365 days is {final_results['PFOA']:.2f} ng/day.")


if __name__ == '__main__':
    calculate_chemical_accumulation()