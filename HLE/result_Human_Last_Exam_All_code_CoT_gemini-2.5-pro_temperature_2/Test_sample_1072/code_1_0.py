# --- Constants and Parameters ---
T_surface = 25.0  # Surface temperature in °C
geo_grad = 25.0  # Geothermal gradient in °C/km
T_ZHe = 180.0  # Closure temperature for Zircon (U-Th)/He in °C
T_AHe = 70.0   # Closure temperature for Apatite (U-Th)/He in °C
total_time = 100.0 # Total duration of exhumation in Myr

print("--- Age Calculation for Statement [H] ---")

# --- Sample 1: Zircon from Pluton ---
print("\n[Sample 1: Zircon]")
depth_1 = 15.0 # km
T1_start = T_surface + depth_1 * geo_grad
rate1 = (T1_start - T_surface) / total_time
time_to_cool1 = (T1_start - T_ZHe) / rate1
date1 = total_time - time_to_cool1

print(f"Starting temperature at 100 Ma = {T_surface}°C + {depth_1} km * {geo_grad}°C/km = {T1_start}°C")
print(f"Cooling rate = ({T1_start} - {T_surface})°C / {total_time} Myr = {rate1:.2f} °C/Myr")
print(f"Time to cool to ZHe closure temp ({T_ZHe}°C) = ({T1_start} - {T_ZHe}) / {rate1:.2f} = {time_to_cool1:.2f} Myr")
print(f"Final Date = {total_time} Ma - {time_to_cool1:.2f} Myr = {date1:.2f} Ma")

# --- Sample 2: Apatite from Sediments ---
print("\n[Sample 2: Apatite]")
T2_start = 250.0 # °C
rate2 = (T2_start - T_surface) / total_time
time_to_cool2 = (T2_start - T_AHe) / rate2
date2 = total_time - time_to_cool2

print(f"Starting temperature at 100 Ma = {T2_start}°C")
print(f"Cooling rate = ({T2_start} - {T_surface})°C / {total_time} Myr = {rate2:.2f} °C/Myr")
print(f"Time to cool to AHe closure temp ({T_AHe}°C) = ({T2_start} - {T_AHe}) / {rate2:.2f} = {time_to_cool2:.2f} Myr")
print(f"Final Date = {total_time} Ma - {time_to_cool2:.2f} Myr = {date2:.2f} Ma")


# --- Sample 3: Apatite from Rhyolite ---
print("\n[Sample 3: Rhyolite]")
date3 = 90.0 # Ma
print(f"Erupted and cooled rapidly at 90 Ma.")
print(f"Final Date = {date3:.2f} Ma")

# --- Comparison ---
print("\n--- Age Comparison ---")
print(f"Ages: Sample 1 ≈ {date1:.0f} Ma, Sample 2 ≈ {date2:.0f} Ma, Sample 3 = {date3:.0f} Ma")
if date3 > date1 and date1 > date2:
    print("Conclusion: Sample 3 is oldest, Sample 2 is youngest. Statement [H] is TRUE.")
else:
    print("Conclusion: Age ordering does not match statement [H].")
