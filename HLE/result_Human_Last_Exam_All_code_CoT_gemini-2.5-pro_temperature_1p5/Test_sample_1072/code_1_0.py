import sys
import io

# Redirect stdout to capture output
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

# --- Step 1: Define Constants and Parameters ---

# Physical constants
T_surface = 25.0  # Surface temperature in °C
geo_gradient = 25.0  # Geothermal gradient in °C/km
Tc_zircon = 190.0  # Approximate closure temperature for Zircon (U-Th)/He in °C
Tc_apatite = 75.0  # Approximate closure temperature for Apatite (U-Th)/He in °C

# Sample parameters
start_time = 100.0 # Ma

# --- Step 2: Calculate Age of Sample 1 (Zircon) ---
print("--- Analysis of Sample 1 (Zircon) ---")
s1_start_depth = 15.0  # km
s1_start_time = 100.0  # Ma
s1_start_temp = T_surface + s1_start_depth * geo_gradient
print(f"Initial temperature at {s1_start_time} Ma: {T_surface}°C + {s1_start_depth} km * {geo_gradient}°C/km = {s1_start_temp}°C (well above Zircon Tc of {Tc_zircon}°C)")

s1_exhumation_rate = s1_start_depth / s1_start_time  # km/Myr
print(f"Exhumation rate: {s1_start_depth} km / {s1_start_time} Myr = {s1_exhumation_rate:.2f} km/Myr")

s1_closure_depth = (Tc_zircon - T_surface) / geo_gradient
print(f"Depth of Zircon closure ({Tc_zircon}°C): ({Tc_zircon}°C - {T_surface}°C) / {geo_gradient}°C/km = {s1_closure_depth:.2f} km")

s1_exhumed_dist_to_closure = s1_start_depth - s1_closure_depth
time_to_cool = s1_exhumed_dist_to_closure / s1_exhumation_rate

s1_age = s1_start_time - time_to_cool
print(f"Equation for Age 1: {s1_start_time} Ma - (({s1_start_depth} km - {s1_closure_depth:.2f} km) / {s1_exhumation_rate:.2f} km/Myr)")
print(f"Calculated Age for Sample 1: {s1_age:.1f} Ma\n")


# --- Step 3: Calculate Age of Sample 2 (Apatite) ---
print("--- Analysis of Sample 2 (Apatite) ---")
s2_start_time = 100.0 # Ma
s2_start_temp = 250.0 # °C
print(f"Thermal clock reset at {s2_start_time} Ma to {s2_start_temp}°C (well above Apatite Tc of {Tc_apatite}°C)")

s2_start_depth = (s2_start_temp - T_surface) / geo_gradient
print(f"Initial depth at {s2_start_time} Ma: ({s2_start_temp}°C - {T_surface}°C) / {geo_gradient}°C/km = {s2_start_depth:.2f} km")

s2_exhumation_rate = s2_start_depth / s2_start_time
print(f"Exhumation rate: {s2_start_depth:.2f} km / {s2_start_time} Myr = {s2_exhumation_rate:.2f} km/Myr")

s2_closure_depth = (Tc_apatite - T_surface) / geo_gradient
print(f"Depth of Apatite closure ({Tc_apatite}°C): ({Tc_apatite}°C - {T_surface}°C) / {geo_gradient}°C/km = {s2_closure_depth:.2f} km")

s2_exhumed_dist_to_closure = s2_start_depth - s2_closure_depth
time_to_cool_s2 = s2_exhumed_dist_to_closure / s2_exhumation_rate

s2_age = s2_start_time - time_to_cool_s2
print(f"Equation for Age 2: {s2_start_time} Ma - (({s2_start_depth:.2f} km - {s2_closure_depth:.2f} km) / {s2_exhumation_rate:.2f} km/Myr)")
print(f"Calculated Age for Sample 2: {s2_age:.1f} Ma\n")


# --- Step 4: State Age of Sample 3 (Apatite) ---
print("--- Analysis of Sample 3 (Apatite) ---")
s3_age = 90.0 # Ma
print(f"Sample 3 is a volcanic rhyolite erupted at {s3_age:.1f} Ma.")
print("It cooled almost instantaneously, so its Apatite (U-Th)/He age is the eruption age.")
print(f"Age for Sample 3: {s3_age:.1f} Ma\n")

# --- Step 5: Evaluate Age Statements ---
print("--- Comparing Ages ---")
print(f"Sample 1 Age: ~{s1_age:.0f} Ma")
print(f"Sample 2 Age: ~{s2_age:.0f} Ma")
print(f"Sample 3 Age: ~{s3_age:.0f} Ma")
print(f"The age order is: Sample 3 ({s3_age:.0f} Ma) > Sample 1 ({s1_age:.0f} Ma) > Sample 2 ({s2_age:.0f} Ma).")
print("Therefore, 'Sample 3 dates are oldest and sample 2 dates are youngest.'")
print("This means statement [H] is TRUE.\n")

# --- Step 6: Evaluate Correlation Statements ---
print("--- Evaluating Correlations ---")

# Date-eU Correlations
print("Statement [B] vs [A]: Sample 1 (Zircon) date-eU correlation")
print("Sample 1 cooled slowly. In robust minerals like Zircon, higher eU (radiation damage) creates more traps for He, increasing He retention and the effective closure temperature. This leads to an older age. This is a POSITIVE date-eU correlation. Statement [B] is TRUE.")

print("\nStatement [C] vs [D]: Sample 2 (Apatite) date-eU correlation")
print("Sample 2 also cooled slowly. However, Apatite is more susceptible to radiation damage. In cases of slow cooling and high damage accumulation, the crystal lattice can become so damaged that it forms interconnected fast-diffusion pathways. This would lead to lower He retention and a younger age for high-eU grains. This is a NEGATIVE date-eU correlation. Given the provided answer choices, this nuanced interpretation is the most likely intended logic. Statement [C] is TRUE.")

# Date-Radius Correlations
print("\nStatements [E] vs [F] vs [G]: Date-Radius correlation")
print("For both Samples 1 and 2, which underwent slow cooling, the date-radius correlation is governed by diffusion principles. Larger grains have longer diffusion path lengths, making it harder for He to escape. This leads to better He retention, a higher effective closure temperature, and an older age. Therefore, both samples should exhibit a POSITIVE date-radius correlation. Statement [E] is TRUE.\n")

# --- Step 7: Final Conclusion ---
print("--- Final Conclusion ---")
print("The true statements are:")
print("[B] Sample 1 has a positive date-eU correlation.")
print("[C] Sample 2 has a negative date-eU correlation.")
print("[E] Samples 1 and 2 have a positive date-radius correlation")
print("[H] Sample 3 dates are oldest and sample 2 dates are youngest.")
print("This combination {B, C, E, H} corresponds to answer choice I.")

# Restore stdout and print the captured output
sys.stdout = old_stdout
print(captured_output.getvalue())
print("<<<I>>>")