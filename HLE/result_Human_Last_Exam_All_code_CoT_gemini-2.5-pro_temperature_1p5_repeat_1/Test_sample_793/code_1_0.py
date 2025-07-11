import numpy as np
from scipy.optimize import minimize

def solve_contaminant_transport():
    """
    Solves the contaminant transport problem as described.
    Calculates peak concentrations, identifies the highest one, and assesses mixture toxicity.
    """
    # 1. Parameter Extraction and Calculation
    # Environmental parameters
    x = 100.0  # m, distance to spring
    v = 0.5  # m/d, advection velocity
    alpha = 0.5  # m, dispersivity factor
    D = alpha * v  # m²/d, dispersion coefficient
    
    # Fracture and coal seam properties
    fracture_length = 80.0 # m
    fracture_height = 10.0 # m
    fracture_area = fracture_length * fracture_height # m²
    
    coal_density = 1346.0  # kg/m³, particle density
    water_content = 0.30  # porosity (n)
    foc = 0.50  # fraction of organic carbon
    bulk_density = coal_density * (1 - water_content)  # kg/m³

    # Injected water volume
    total_water_litres = 1_500_000.0

    # Chemical-specific data
    chemicals = {
        'Atrazine': {
            'product_fraction': 0.01,
            'conc_in_product': 40.0,  # μg/L
            'log_koc': 2.20,
            'half_life_days': 90.0,
            'ec50': 100.0  # μg/L
        },
        'PFOS': {
            'product_fraction': 0.001,
            'conc_in_product': 300.0, # μg/L
            'log_koc': 3.65,
            'half_life_days': 14965.0,
            'ec50': 480.0  # μg/L
        },
        'Endosulfan': {
            'product_fraction': 0.005,
            'conc_in_product': 20.0,  # μg/L
            'log_koc': 4.3,
            'half_life_days': 60.0,
            'ec50': 560.0  # μg/L
        }
    }

    # 2. Calculate derived parameters for each chemical
    for name, data in chemicals.items():
        # Mass per unit area (MTotal) in μg/m²
        total_mass_ug = total_water_litres * data['product_fraction'] * data['conc_in_product']
        data['MTotal'] = total_mass_ug / fracture_area
        
        # Decay rate (k) in 1/d
        data['k'] = np.log(2) / data['half_life_days']
        
        # Partitioning coefficient (Kd) in L/kg
        koc = 10**data['log_koc']
        data['Kd'] = koc * foc
        
        # Retardation factor (R), dimensionless
        # R = 1 + (ρ_b / n) * Kd, converting Kd from L/kg to m³/kg
        data['R'] = 1 + (bulk_density / water_content) * (data['Kd'] / 1000.0)

    # Advection-Diffusion Equation C(x,t)
    def calculate_concentration(t, MTotal, k, R):
        # Prevent division by zero or log of zero for t=0
        if t <= 1e-9:
            return 0.0
        
        # All units are consistent (m, d, kg, L, μg)
        # Result C is in μg/m³
        term1 = MTotal / np.sqrt(4 * np.pi * D * t)
        exponent_arg = -((x - v * t)**2) / (4 * D * t)
        term2 = np.exp(exponent_arg)
        term3 = np.exp(-k * t)
        retardation_factor = 1.0 / R
        
        concentration_ug_m3 = term1 * term2 * term3 * retardation_factor
        # Convert concentration from μg/m³ to μg/L
        return concentration_ug_m3 / 1000.0

    print("Calculating Peak Contaminant Concentrations at the Spring (x=100m)")
    print("=" * 65)

    results = {}
    sum_rq = 0.0
    
    # 3. Find peak concentration for each chemical
    for name, data in chemicals.items():
        # Objective function to minimize (negative of concentration)
        objective_func = lambda t: -calculate_concentration(t[0], data['MTotal'], data['k'], data['R'])
        
        # Find time of peak by minimizing the negative of the concentration function
        # Initial guess for time is t = x/v = 100/0.5 = 200 days
        opt_result = minimize(objective_func, x0=[200], bounds=[(1e-6, None)])
        
        peak_time = opt_result.x[0]
        peak_conc = -opt_result.fun # convert back to positive
        
        results[name] = {
            'peak_conc': peak_conc,
            'peak_time': peak_time
        }
        
        # Output the calculation with numbers plugged in
        print(f"Calculation for {name}:")
        print(f"Equation: C(x,t) = MTotal/sqrt(4πDt) * exp(-(x-vt)²/4Dt) * exp(-kt) * 1/R")
        print(f"At peak time t = {peak_time:.2f} days:")
        
        # Formatted string showing the equation with substituted values
        equation_str = (
            f"  C({x}, {peak_time:.2f}) = ({data['MTotal']:.2e} μg/m² / sqrt(4π * {D:.2f} * {peak_time:.2f})) * "
            f"exp(-({x} - {v}*{peak_time:.2f})² / (4*{D:.2f}*{peak_time:.2f})) * "
            f"exp(-{data['k']:.5f}*{peak_time:.2f}) * (1/{data['R']:.2f})"
        )
        print(equation_str)
        print(f"  Resulting Peak Concentration = {peak_conc:.6f} µg/L\n")
        
        # 4. Calculate Risk Quotient and add to sum
        rq = peak_conc / data['ec50']
        sum_rq += rq

    # Determine highest concentration
    highest_chem = max(results, key=lambda name: results[name]['peak_conc'])
    max_conc = results[highest_chem]['peak_conc']
    
    # Assess mixture effect based on sum of RQs
    # Standard assumption is additivity.
    mixture_effect_type = "additive"
    
    print("=" * 65)
    print("Summary of Results:")
    print(f"\nThe highest concentration of an individual chemical reaching the spring is from {highest_chem}.")
    print(f"Highest Concentration: {max_conc:.6f} µg/L")
    print(f"\nMixture Effect Analysis (assuming additivity):")
    print(f"Sum of Risk Quotients (C_peak / EC50) = {sum_rq:.6f}")
    
    # Final answer in requested format
    print(f"<<<{max_conc:.6f}, {mixture_effect_type}>>>")

if __name__ == '__main__':
    solve_contaminant_transport()