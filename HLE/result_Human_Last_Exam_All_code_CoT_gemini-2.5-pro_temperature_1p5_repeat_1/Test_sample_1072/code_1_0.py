import sys
import io

# Capture original stdout
original_stdout = sys.stdout
# Create a string buffer
string_buffer = io.StringIO()
# Redirect stdout to the string buffer
sys.stdout = string_buffer

# --- Problem Setup ---
# This script calculates the expected (U-Th)/He ages for three different geological scenarios
# to determine their relative order and evaluate statement [H]. It also includes the reasoning
# for the expected date-eU and date-radius correlations to select the final correct answer choice.

# --- Part 1: Age Calculation ---

# Define geological constants
SURFACE_T = 25      # Surface temperature in °C
GEO_GRADIENT = 25   # Geothermal gradient in °C/km
PRESENT_DAY = 0     # Age of the present in Ma

# Define mineral and sample-specific parameters
# Closure Temperature (Tc) is the approximate temperature below which a mineral
# effectively retains Helium-4.
S1 = {'name': 'Sample 1 (Zircon Pluton)', 'Tc': 180, 'start_time': 100, 'start_depth': 15}
S2 = {'name': 'Sample 2 (Apatite Sediment)', 'Tc': 70, 'start_time': 100, 'start_temp': 250}
S3 = {'name': 'Sample 3 (Apatite Rhyolite)', 'age': 90}

print("### Step 1: Calculating the Expected Ages ###")

# --- Calculation for Sample 1 ---
print(f"\n--- Analysis for {S1['name']} ---")
# The sample underwent steady exhumation for 100 Myr from a depth of 15 km.
exhumation_duration_s1 = S1['start_time'] - PRESENT_DAY
exhumation_rate_s1 = S1['start_depth'] / exhumation_duration_s1
print(f"The exhumation rate is {S1['start_depth']} km / {exhumation_duration_s1} Myr = {exhumation_rate_s1:.2f} km/Myr.")
# The He clock starts when the sample cools below its closure temperature (Tc).
# First, find the depth corresponding to Tc.
# Equation: T(depth) = T_surface + gradient * depth
# Rearranged: depth = (T(depth) - T_surface) / gradient
closure_depth_s1 = (S1['Tc'] - SURFACE_T) / GEO_GRADIENT
print(f"The closure depth for Zircon (Tc={S1['Tc']}°C) is ({S1['Tc']} - {SURFACE_T}) / {GEO_GRADIENT} = {closure_depth_s1:.2f} km.")
# Next, calculate how long it took to travel from the start depth to the closure depth.
distance_exhumed_s1 = S1['start_depth'] - closure_depth_s1
time_to_closure_s1 = distance_exhumed_s1 / exhumation_rate_s1
print(f"The time to exhume to this depth is ({S1['start_depth']} - {closure_depth_s1:.2f}) km / {exhumation_rate_s1:.2f} km/Myr = {time_to_closure_s1:.2f} Myr.")
# The age is the start time minus the time elapsed until closure.
age_s1 = S1['start_time'] - time_to_closure_s1
print(f"The final predicted age is {S1['start_time']} Ma - {time_to_closure_s1:.2f} Ma = {age_s1:.2f} Ma.")

# --- Calculation for Sample 2 ---
print(f"\n--- Analysis for {S2['name']} ---")
# The sample was heated to 250°C at 100 Ma, which reset its He clock. It then exhumed.
# First, find the starting depth at 100 Ma from the temperature.
start_depth_s2 = (S2['start_temp'] - SURFACE_T) / GEO_GRADIENT
print(f"The initial depth at 100 Ma (T={S2['start_temp']}°C) was ({S2['start_temp']} - {SURFACE_T}) / {GEO_GRADIENT} = {start_depth_s2:.1f} km.")
# Now, calculate its exhumation rate.
exhumation_duration_s2 = S2['start_time'] - PRESENT_DAY
exhumation_rate_s2 = start_depth_s2 / exhumation_duration_s2
print(f"The exhumation rate is {start_depth_s2:.1f} km / {exhumation_duration_s2} Myr = {exhumation_rate_s2:.2f} km/Myr.")
# Find the closure depth for Apatite.
closure_depth_s2 = (S2['Tc'] - SURFACE_T) / GEO_GRADIENT
print(f"The closure depth for Apatite (Tc={S2['Tc']}°C) is ({S2['Tc']} - {SURFACE_T}) / {GEO_GRADIENT} = {closure_depth_s2:.2f} km.")
# Calculate the time to exhume to this closure depth.
distance_exhumed_s2 = start_depth_s2 - closure_depth_s2
time_to_closure_s2 = distance_exhumed_s2 / exhumation_rate_s2
print(f"The time to exhume to this depth is ({start_depth_s2:.1f} - {closure_depth_s2:.2f}) km / {exhumation_rate_s2:.2f} km/Myr = {time_to_closure_s2:.2f} Myr.")
# Calculate the final age.
age_s2 = S2['start_time'] - time_to_closure_s2
print(f"The final predicted age is {S2['start_time']} Ma - {time_to_closure_s2:.2f} Ma = {age_s2:.2f} Ma.")

# --- Age for Sample 3 ---
print(f"\n--- Analysis for {S3['name']} ---")
age_s3 = S3['age']
print(f"This sample erupted at {age_s3} Ma. Volcanic eruption involves rapid cooling, so the He age is the eruption age.")
print(f"The final predicted age is {age_s3:.2f} Ma.")

# --- Conclusion on Ages ---
print("\n### Conclusion on Relative Ages ###")
print(f"Sample 1 Age ≈ {age_s1:.0f} Ma")
print(f"Sample 2 Age ≈ {age_s2:.0f} Ma")
print(f"Sample 3 Age ≈ {age_s3:.0f} Ma")
print(f"Comparing the ages: {age_s3:.0f} Ma (Sample 3) > {age_s1:.0f} Ma (Sample 1) > {age_s2:.0f} Ma (Sample 2).")
print("This means: Sample 3 dates are oldest and sample 2 dates are youngest.")
print("Therefore, statement [H] is TRUE.")

# --- Part 2: Correlation Analysis and Final Answer Selection ---
print("\n### Step 2: Analyzing Correlations and Selecting the Final Answer ###")
print("\n[E], [F], [G] Date-Radius Correlation:")
print(" - For slowly cooled samples (1 and 2), larger crystals retain He better and thus have higher effective closure temperatures, yielding older ages. This predicts a POSITIVE date-radius correlation for both.")
print(" - This means statement [E] (Samples 1 and 2 have a positive date-radius correlation) is TRUE.")

print("\n[A], [B], [C], [D] Date-eU Correlation:")
print(" - For slowly cooled samples, higher eU (Uranium) causes more radiation damage, which inhibits He diffusion, leading to older ages. This predicts a POSITIVE date-eU correlation.")
print(" - [B] Sample 1 (Zircon) has a positive date-eU correlation. This is the standard, robust model for zircon. TRUE.")
print(" - [D] Sample 2 (Apatite) has a positive date-eU correlation. This is the standard model for apatite. TRUE.")

print("\n### Synthesis ###")
print("Based on standard models, the true statements are [B], [D], [E], and [H].")
print("This combination {B, D, E, H} is not an available answer choice.")
print("This implies we must choose the 'best-fit' answer, which may rely on a non-standard exception.")
print("Let's re-examine the choices that contain our most certain conclusions: [H] and [E]. These are choices H and I.")
print(" - Choice H: {A, D, E, H}. This would mean Zircon has a negative date-eU correlation [A]. This contradicts the standard model.")
print(" - Choice I: {B, C, E, H}. This would mean Apatite has a negative date-eU correlation [C].")
print("Between the two, it is more plausible for the complex Apatite system (with damage annealing) to exhibit exceptional behavior than the more robust Zircon system.")
print("Therefore, the most logical choice is to assume the question implies an exception for Apatite, making [C] true in this context.")
print("This makes the final set of correct statements {B, C, E, H}.")

# Restore original stdout
sys.stdout = original_stdout
# Get the captured output
output = string_buffer.getvalue()

# Print the final output
print(output)
print("The best fit answer is I, which corresponds to statements B, C, E, and H being true.")

<<<I>>>