# Patient-specific variables
bsa_m2 = 0.95  # Body Surface Area in square meters

# Standard constants for the BSA method
daily_fluid_per_bsa = 1500  # Standard maintenance fluid rate in mL per m^2 per day
hours_in_a_day = 24

# --- Calculation ---

# 1. Calculate the total daily maintenance fluid requirement in mL
total_daily_fluid = daily_fluid_per_bsa * bsa_m2

# 2. Convert the daily requirement to an hourly rate in cc/hr (1 mL = 1 cc)
hourly_rate = total_daily_fluid / hours_in_a_day

# --- Output the results ---
print("Calculating Pediatric Maintenance Fluids using the BSA Method")
print(f"Formula: (1500 mL/m²/day * BSA) / 24 hours")
print("-" * 30)
print("Step 1: Calculate total daily fluid volume.")
print(f"{daily_fluid_per_bsa} mL/m²/day * {bsa_m2} m² = {total_daily_fluid:.2f} mL/day")
print("\nStep 2: Convert daily volume to hourly rate.")
print(f"{total_daily_fluid:.2f} mL/day / {hours_in_a_day} hours = {hourly_rate:.2f} cc/hr")
print("\nFinal Answer:")
print(f"The patient's maintenance fluid rate should be {hourly_rate:.2f} cc/hr.")
print(f"The full calculation is: ({daily_fluid_per_bsa} * {bsa_m2}) / {hours_in_a_day} = {hourly_rate:.2f}")

final_answer = round(hourly_rate, 2)