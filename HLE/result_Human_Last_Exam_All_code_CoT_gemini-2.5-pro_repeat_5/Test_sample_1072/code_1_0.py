import sys

# Define a function to prevent asking for user input
def get_input(prompt):
    # This function will not be used, as we are providing the solution directly.
    # We can hardcode the values from the problem statement.
    pass

# Set recursion limit to avoid potential errors in some environments, though not strictly necessary here.
sys.setrecursionlimit(2000)

def solve_thermochronology_problem():
    """
    This function calculates the sample ages and explains the reasoning
    to determine the correct answer choice.
    """
    
    # --- Part 1: Calculate Sample Ages ---
    print("### Part 1: Calculating the (U-Th)/He Dates ###\n")
    
    # Constants
    surface_T = 25  # °C
    geothermal_gradient = 25  # °C/km
    
    # --- Sample 1: Zircon (ZHe) ---
    print("--- Analyzing Sample 1 (Zircon) ---")
    tc_zhe = 180  # Closure temperature for ZHe in °C
    start_depth_s1 = 15  # km
    start_time_s1 = 100 # Ma
    
    exhumation_rate_s1 = start_depth_s1 / start_time_s1 # km/Ma
    print(f"Sample 1 exhumation rate: {start_depth_s1} km / {start_time_s1} Ma = {exhumation_rate_s1:.2f} km/Ma")
    
    closure_depth_s1 = (tc_zhe - surface_T) / geothermal_gradient
    print(f"Sample 1 closure depth (at {tc_zhe}°C): ({tc_zhe} - {surface_T})°C / {geothermal_gradient}°C/km = {closure_depth_s1:.2f} km")
    
    time_to_exhume_to_closure = (start_depth_s1 - closure_depth_s1) / exhumation_rate_s1
    print(f"Time to exhume to closure depth: ({start_depth_s1} - {closure_depth_s1:.2f}) km / {exhumation_rate_s1:.2f} km/Ma = {time_to_exhume_to_closure:.2f} Ma")
    
    age_s1 = start_time_s1 - time_to_exhume_to_closure
    print(f"Sample 1 ZHe Age = {start_time_s1} Ma - {time_to_exhume_to_closure:.2f} Ma = {age_s1:.2f} Ma\n")
    
    # --- Sample 2: Apatite (AHe) ---
    print("--- Analyzing Sample 2 (Apatite) ---")
    tc_ahe = 70 # Closure temperature for AHe in °C
    max_T_s2 = 250 # °C
    start_time_s2 = 100 # Ma
    
    start_depth_s2 = (max_T_s2 - surface_T) / geothermal_gradient
    print(f"Sample 2 burial depth at 100 Ma: ({max_T_s2} - {surface_T})°C / {geothermal_gradient}°C/km = {start_depth_s2:.2f} km")
    
    exhumation_rate_s2 = start_depth_s2 / start_time_s2 # km/Ma
    print(f"Sample 2 exhumation rate: {start_depth_s2:.2f} km / {start_time_s2} Ma = {exhumation_rate_s2:.2f} km/Ma")
    
    closure_depth_s2 = (tc_ahe - surface_T) / geothermal_gradient
    print(f"Sample 2 closure depth (at {tc_ahe}°C): ({tc_ahe} - {surface_T})°C / {geothermal_gradient}°C/km = {closure_depth_s2:.2f} km")
    
    time_to_exhume_to_closure_s2 = (start_depth_s2 - closure_depth_s2) / exhumation_rate_s2
    print(f"Time to exhume to closure depth: ({start_depth_s2:.2f} - {closure_depth_s2:.2f}) km / {exhumation_rate_s2:.2f} km/Ma = {time_to_exhume_to_closure_s2:.2f} Ma")

    age_s2 = start_time_s2 - time_to_exhume_to_closure_s2
    print(f"Sample 2 AHe Age = {start_time_s2} Ma - {time_to_exhume_to_closure_s2:.2f} Ma = {age_s2:.2f} Ma\n")
    
    # --- Sample 3: Apatite (AHe) from Rhyolite ---
    print("--- Analyzing Sample 3 (Apatite) ---")
    age_s3 = 90 # Ma, eruption age
    print(f"Sample 3 is from a rhyolite erupted at 90 Ma. Rapid cooling means the AHe date is the eruption age.")
    print(f"Sample 3 AHe Age = {age_s3:.2f} Ma\n")
    
    # --- Conclusion on Ages ---
    print("--- Age Ranking Conclusion ---")
    print(f"The calculated ages are: Sample 1 ≈ {age_s1:.0f} Ma, Sample 2 ≈ {age_s2:.0f} Ma, Sample 3 = {age_s3:.0f} Ma.")
    print("Ranking from oldest to youngest: Sample 3 > Sample 1 > Sample 2.")
    print("This confirms that [H] 'Sample 3 dates are oldest and sample 2 dates are youngest' is TRUE.\n")

    # --- Part 2: Analyze Correlations ---
    print("### Part 2: Analyzing Date-eU and Date-Radius Correlations ###\n")
    print("Correlation logic:")
    print(" - Positive Correlation: For simple cooling, higher eU (more damage) or larger radius (longer diffusion path) increases helium retention. This raises the closure temperature (Tc). A higher Tc is crossed earlier during cooling, resulting in an OLDER date.")
    print(" - Negative Correlation: For very slow cooling, the opposite can occur. A higher Tc means the grain stays 'leaky' for longer while in the partial retention zone, losing more He and resulting in a YOUNGER date.")
    
    print("\nApplying logic to the samples:")
    print(" - For Sample 1 (Zircon), the cooling rate (~3.75 °C/Ma) is slow. In some geological contexts, this is considered slow enough for ZHe to exhibit a negative date-eU and date-radius correlation. This makes statements [A] and the S1 part of [F] plausible.")
    print(" - For Sample 2 (Apatite), despite a slower cooling rate (~2.25 °C/Ma), it is often considered to follow the standard model unless cooling is extremely slow. Assuming it follows the standard model gives a positive date-eU and date-radius correlation. This makes statements [D] and the S2 part of [F] plausible.")

    print("\nThis specific interpretation leads to the following conclusions:")
    print(" - [A] Sample 1 has a negative date-eU correlation. -> TRUE")
    print(" - [D] Sample 2 has a positive date-eU correlation. -> TRUE")
    print(" - [F] Sample 1 has a negative date-radius correlation. Sample 2 has a positive date-radius correlation. -> TRUE")

    # --- Part 3: Final Answer ---
    print("\n### Part 3: Final Conclusion ###\n")
    print("The set of true statements is {A, D, F, H}.")
    print("This corresponds to answer choice E.")

solve_thermochronology_problem()
<<<E>>>