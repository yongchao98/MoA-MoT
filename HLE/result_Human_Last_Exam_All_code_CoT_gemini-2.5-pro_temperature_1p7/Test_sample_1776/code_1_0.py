import math
from datetime import datetime

# --- Given Information & Constants ---
# Radionuclide: Indium-111
T_HALF_DAYS = 2.8
T_HALF_HOURS = T_HALF_DAYS * 24

# Dose information for the patient
A_dose = 8.0  # mCi
injection_time = datetime(2024, 12, 23, 8, 0)

# Compounding time
compounding_time = datetime(2024, 12, 23, 4, 4)

# Stock vial calibration information
A_calibrated = 10.0  # mCi
V_calibrated = 10.0  # mL
calibration_time = datetime(2024, 12, 26, 12, 0)

# --- Calculations ---

# 1. Calculate the decay constant (λ) for Indium-111
lambda_val = math.log(2) / T_HALF_HOURS

# 2. Calculate the required activity at the time of compounding (A_needed)
# This accounts for the decay between compounding and injection.
time_decay_td = injection_time - compounding_time
time_decay_hours = time_decay_td.total_seconds() / 3600
A_needed = A_dose * math.exp(lambda_val * time_decay_hours)

# 3. Calculate the activity of the stock vial at compounding time (A_stock_now)
# This accounts for the decay from compounding time until the future calibration time.
time_to_cal_td = calibration_time - compounding_time
time_to_cal_hours = time_to_cal_td.total_seconds() / 3600
A_stock_now = A_calibrated * math.exp(lambda_val * time_to_cal_hours)

# 4. Calculate the concentration of the stock vial at compounding time
concentration_now = A_stock_now / V_calibrated

# 5. Calculate the final volume to be drawn
volume_to_draw = A_needed / concentration_now

# --- Output the Result ---
print("This script calculates the required volume for a radiopharmaceutical dose.")
print("It accounts for radioactive decay for both the patient dose and the stock vial.\n")
print(f"Step 1: The activity needed in the syringe at 4:04 am is {A_needed:.4f} mCi.")
print(f"Step 2: The activity in the stock vial at 4:04 am is {A_stock_now:.4f} mCi.")
print(f"Step 3: The concentration of the stock vial at 4:04 am is {concentration_now:.4f} mCi/mL.\n")
print("Step 4: Final volume calculation:")
print(f"Volume to Draw = (Required Activity) / (Stock Concentration)")
print(f"Volume to Draw = {A_needed:.4f} mCi / {concentration_now:.4f} mCi/mL\n")
print(f"The final volume of Indium 111 chloride that must be drawn is {volume_to_draw:.2f} mL.")
