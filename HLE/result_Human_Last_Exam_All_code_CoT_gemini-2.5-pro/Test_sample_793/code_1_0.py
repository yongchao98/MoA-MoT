import numpy as np
from scipy.optimize import minimize_scalar

def solve_transport_problem():
    """
    Calculates the highest contaminant concentration reaching a spring and assesses the mixture effect.
    """
    # Step 1: Define general parameters
    x = 100.0  # distance to spring [m]
    v = 0.5  # advection velocity [m/d]
    dispersivity = 0.5  # dispersivity [m]
    D = v * dispersivity  # dispersion coefficient [m^2/d]
    rho_b = 1346.0  # coal seam bulk density [kg/m^3]
    n = 0.30  # porosity (water content) [dimensionless]
    f_oc = 0.50  # fraction of organic carbon [dimensionless]
    A = 80.0 * 10.0 # cross-sectional area of flow, from fracture dimensions [m^2]
    injected_volume_L = 1500000.0  # [L]

    # Step 2: Define chemical-specific data
    chemicals = [
        {
            "name": "Atrazine",
            "product_conc_ug_L": 40, "product_fraction": 0.01,
            "log_koc": 2.20, "half_life": 90, "ec50": 100
        },
        {
            "name": "PFOS",
            "product_conc_ug_L": 300, "product_fraction": 0.001,
            "log_koc": 3.65, "half_life": 14965, "ec50": 480
        },
        {
            "name": "Endosulfan",
            "product_conc_ug_L": 20, "product_fraction": 0.005,
            "log_koc": 4.3, "half_life": 60, "ec50": 560
        }
    ]

    # Calculate derived parameters for each chemical
    for chem in chemicals:
        conc_in_water_ug_L = chem["product_conc_ug_L"] * chem["product_fraction"]
        chem["M"] = conc_in_water_ug_L * injected_volume_L  # Total mass [ug]
        chem["k"] = np.log(2) / chem["half_life"]  # Decay rate [1/d]
        koc = 10**chem["log_koc"]
        kd_L_kg = koc * f_oc
        kd_m3_kg = kd_L_kg / 1000.0  # Convert L/kg to m^3/kg
        chem["R"] = 1 + (rho_b / n) * kd_m3_kg  # Retardation factor

    # Step 3: Define concentration function
    def calculate_C(t, chem):
        if t <= 0: return 0.0
        # Equation: C = (M / (A*n*R*sqrt(4*pi*D*t))) * exp(-((x-v*t)**2)/(4*D*t)) * exp(-k*t)
        M_term = chem["M"] / (A * n * chem["R"])
        sqrt_term = np.sqrt(4 * np.pi * D * t)
        gauss_term = np.exp(-((x - v * t)**2) / (4 * D * t))
        decay_term = np.exp(-chem["k"] * t)
        
        conc_ug_m3 = (M_term / sqrt_term) * gauss_term * decay_term
        return conc_ug_m3

    # Step 4: Find the peak concentration for each chemical
    results = []
    print("--- Finding Peak Concentrations for Each Chemical ---\n")
    for chem in chemicals:
        # Objective function to minimize (negative concentration)
        def neg_C(t):
            return -calculate_C(t, chem)

        # Find the time of peak concentration by minimizing the negative concentration
        res = minimize_scalar(neg_C, bounds=(1, 1000), method='bounded')
        
        t_peak = res.x
        c_peak_ug_m3 = -res.fun
        c_peak_ug_L = c_peak_ug_m3 / 1000.0 # Convert ug/m^3 to ug/L
        
        results.append({
            "name": chem["name"], "c_peak_ug_L": c_peak_ug_L, "t_peak_days": t_peak
        })
        print(f"{chem['name']}:")
        print(f"  - Peak concentration of {c_peak_ug_L:.6f} µg/L occurs at {t_peak:.2f} days.\n")

    # Step 5: Identify the highest concentration and the time it occurs
    highest_res = max(results, key=lambda r: r['c_peak_ug_L'])
    t_at_highest_conc = highest_res['t_peak_days']

    print("--- Highest Concentration Analysis ---\n")
    print(f"The highest individual chemical concentration reaching the spring is from {highest_res['name']}.")
    print(f"Highest Concentration: {highest_res['c_peak_ug_L']:.6f} µg/L")
    print(f"This peak occurs at time t = {t_at_highest_conc:.2f} days.\n")

    # Step 6: Assess the mixture effect at that specific time
    print("--- Mixture Effect Analysis ---\n")
    print(f"To assess the mixture effect, we calculate the concentrations of all chemicals at t = {t_at_highest_conc:.2f} days.")
    
    concentrations_at_peak_time = {}
    for chem in chemicals:
        conc_ug_m3 = calculate_C(t_at_highest_conc, chem)
        conc_ug_L = conc_ug_m3 / 1000.0
        concentrations_at_peak_time[chem['name']] = conc_ug_L
        print(f"  - Concentration of {chem['name']}: {conc_ug_L:.6f} µg/L")
    
    print("\nThe combined effect is assessed using the Hazard Quotient (HQ) method, which assumes concentration additivity.")
    print("The formula is: HQ_mixture = sum(Concentration_i / EC50_i)")
    
    hq_mix = 0
    hq_calcs = []
    for chem in chemicals:
        C_i = concentrations_at_peak_time[chem['name']]
        EC50_i = chem['ec50']
        hq_i = C_i / EC50_i
        hq_mix += hq_i
        hq_calcs.append(f"({C_i:.6f} / {EC50_i})")

    print(f"\nHQ_mixture = {' + '.join(hq_calcs)}")
    print(f"HQ_mixture = {hq_mix:.6f}")
    
    print("\nConclusion on Mixture Effect:")
    print("The Hazard Quotient approach intrinsically assumes an **additive** effect, where the combined toxicity is the sum of the individual toxicities scaled by their potency (EC50). Without specific data on chemical interactions, additivity is the standard assumption.")
    if hq_mix >= 1:
        print("Since the mixture HQ is >= 1, there is a potential risk to the algal community from the combined exposure.")
    else:
        print("Since the mixture HQ is < 1, the risk to the algal community from the combined exposure is considered low.")

    return highest_res['c_peak_ug_L']

if __name__ == '__main__':
    highest_concentration = solve_transport_problem()
    print(f"\n<<<Final Answer Marker>>>")
    print(f"<<<{highest_concentration:.6f}>>>")
