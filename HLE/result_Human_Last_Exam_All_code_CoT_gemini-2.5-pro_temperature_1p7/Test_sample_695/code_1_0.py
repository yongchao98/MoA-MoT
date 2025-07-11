import math

# An explanation of the approach is provided above.
# This script executes the plan to calculate the accumulated mass of PFOS and PFOA in fish.

# --- Given Parameters ---

# Environment
V_water = 10000.0  # L
Qin = 900.0       # L/d
Qout = 1600.0     # L/d
foc = 0.001       # 0.1% as a fraction
t_period = 365.0         # days

# Fish
M_fish = 1000.0   # g
C_food = 100.0    # ng/g
IR_food = 20.0    # g/day
AF_gills = 0.8    # absorption fraction for gills
AF_food = 0.9     # absorption fraction for food
Q_gills = 100.0   # L/day
C_fish_initial = 10.0 # ng/g

# PFOS specific parameters
Cin_pfos = 2.6     # ng/L
t_half_pfos = 91.0 # years
logKow_pfos = 4.0
kelim_pfos = 0.069 # d^-1

# PFOA specific parameters
Cin_pfoa = 211300.0 # ng/L
t_half_pfoa = 238.0 # years
logKow_pfoa = 4.5
kelim_pfoa = 0.023 # d^-1

def solve_accumulation():
    """
    Solves the entire problem by calculating water concentration and then fish accumulation.
    """

    # --- Step 1: Calculate Water Concentration C(t) for each chemical ---

    # The hydraulic residence time rate constant (r) is common for both chemicals
    r = Qout / V_water

    def calculate_steady_state_water_conc(Cin, t_half_years, logKow):
        """
        Calculates the steady-state concentration in the water body using the given formula.
        C(t) = (Cin * Qin / (Qout + Qin * (1 + Kd * foc))) * (1 - e^-(k + r)t)
        As t becomes large, (1 - e^...) approaches 1, giving the steady-state value.
        """
        # Degradation rate constant (k)
        t_half_days = t_half_years * 365.0
        k = math.log(2) / t_half_days
        
        # Partition coefficient between organic carbon and water (Koc)
        logKoc = 0.81 * logKow + 0.01
        Koc = 10**logKoc
        
        # Soil-water partition coefficient (Kd)
        Kd = Koc * foc # The equation seems to use Kd*foc, which might imply Kd=Koc here.
        
        # C_steady_state = (Cin * Qin) / (Qout + Qin * (1 + Koc * foc * foc)) <-- Error in logic
        # The equation states (1 + Kd * foc) where Kd = Koc * foc, so it should be (1 + (Koc*foc)*foc) 
        # A more standard interpretation is that Kd = Koc, which is dimensionally inconsistent.
        # Let's assume the provided expression is literally Kd = a variable to be found from Koc.
        # Let's use Kd = Koc as it seems implied by logKoc = f(logKow).
        # C_steady_state = (Cin * Qin) / (Qout + Qin * (1 + Kd * foc))
        
        numerator = Cin * Qin
        denominator = Qout + Qin * (1 + Koc * foc) # Most likely interpretation of the formula terms
        
        C_ss = numerator / denominator
        return C_ss

    # Calculate steady-state water concentrations for PFOS and PFOA
    C_water_pfos = calculate_steady_state_water_conc(Cin_pfos, t_half_pfos, logKow_pfos)
    C_water_pfoa = calculate_steady_state_water_conc(Cin_pfoa, t_half_pfoa, logKow_pfoa)

    # --- Step 2: Calculate Accumulated Mass in Fish over 365 days ---

    def calculate_accumulated_mass(C_water, kelim):
        """
        Calculates the final concentration in fish and the total accumulated mass.
        """
        # Total uptake rate (U_total) from gills and food, in ng/day
        uptake_gills = C_water * Q_gills * AF_gills
        uptake_food = C_food * IR_food * AF_food
        uptake_total = uptake_gills + uptake_food
        
        # Fish steady-state concentration (C_fish_ss) in ng/g
        # This is the theoretical concentration if exposure is indefinite.
        C_fish_ss = uptake_total / (kelim * M_fish)
        
        # Final fish concentration at t=365 days using the analytical solution:
        # C_fish(t) = C_fish_ss + (C_fish_initial - C_fish_ss) * exp(-kelim * t)
        C_fish_final = C_fish_ss + (C_fish_initial - C_fish_ss) * math.exp(-kelim * t_period)
        
        # Total accumulated mass is the change in total chemical mass in the fish
        accumulated_mass = (C_fish_final - C_fish_initial) * M_fish
        
        return accumulated_mass, C_fish_final

    # Calculate final results for PFOS and PFOA
    accumulated_mass_pfos, C_final_pfos = calculate_accumulated_mass(C_water_pfos, kelim_pfos)
    accumulated_mass_pfoa, C_final_pfoa = calculate_accumulated_mass(C_water_pfoa, kelim_pfoa)

    # --- Step 3: Output the results in the requested format ---

    # PFOS Results
    print("--- PFOS Accumulation Analysis ---")
    print(f"Calculated final fish concentration after {int(t_period)} days (C_final): {C_final_pfos:.4f} ng/g")
    print("\nFinal Equation for PFOS Accumulated Mass:")
    print("Accumulated Mass = (C_final - C_initial) * M_fish")
    print(f"Accumulated Mass = ({C_final_pfos:.4f} ng/g - {C_fish_initial:.4f} ng/g) * {M_fish:.1f} g")
    print(f"Total Accumulated Mass (PFOS) = {accumulated_mass_pfos:.2f} ng\n")

    # PFOA Results
    print("--- PFOA Accumulation Analysis ---")
    print(f"Calculated final fish concentration after {int(t_period)} days (C_final): {C_final_pfoa:.4f} ng/g")
    print("\nFinal Equation for PFOA Accumulated Mass:")
    print("Accumulated Mass = (C_final - C_initial) * M_fish")
    print(f"Accumulated Mass = ({C_final_pfoa:.4f} ng/g - {C_fish_initial:.4f} ng/g) * {M_fish:.1f} g")
    print(f"Total Accumulated Mass (PFOA) = {accumulated_mass_pfoa:.2f} ng")
    
    return accumulated_mass_pfos, accumulated_mass_pfoa

# Run the full calculation and print the results
pfos_mass, pfoa_mass = solve_accumulation()
# The final answer format is specified as <<<answer content>>>.
# As the answer contains two values, I will format them clearly.
# print(f"\n<<<PFOS Accumulated Mass: {pfos_mass:.2f} ng, PFOA Accumulated Mass: {pfoa_mass:.2f} ng>>>")