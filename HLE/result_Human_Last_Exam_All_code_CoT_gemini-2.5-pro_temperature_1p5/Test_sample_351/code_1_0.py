import math

# Patient Information
weight_kg = 25
tbsa_percent = 45
bsa_m2 = 0.95

print("Calculating pediatric maintenance fluid needs for a burn patient.\n")
print(f"Patient data: Weight={weight_kg}kg, TBSA Burn={tbsa_percent}%, BSA={bsa_m2} m^2\n")

# --- Step 1: Calculate Baseline Maintenance Fluids (Holliday-Segar Formula) ---
daily_maintenance_fluid = 0
if weight_kg <= 10:
    daily_maintenance_fluid = weight_kg * 100
elif weight_kg <= 20:
    daily_maintenance_fluid = (10 * 100) + ((weight_kg - 10) * 50)
else:
    daily_maintenance_fluid = (10 * 100) + (10 * 50) + ((weight_kg - 20) * 20)

hourly_maintenance_fluid = daily_maintenance_fluid / 24

print("Step 1: Calculate Baseline Maintenance Rate (Holliday-Segar)")
print(f"Daily baseline maintenance fluid: {daily_maintenance_fluid:.2f} cc/day")
print(f"Hourly baseline maintenance fluid: {daily_maintenance_fluid:.2f} / 24 = {hourly_maintenance_fluid:.2f} cc/hr\n")


# --- Step 2: Calculate Additional Evaporative Fluid Loss ---
burn_area_m2 = bsa_m2 * (tbsa_percent / 100)
# A common constant for evaporative loss is 2500 cc/m^2/day
daily_evaporative_loss = burn_area_m2 * 2500
hourly_evaporative_loss = daily_evaporative_loss / 24

print("Step 2: Calculate Additional Evaporative Loss Rate")
print(f"Burn surface area: {bsa_m2} m^2 * {tbsa_percent / 100} = {burn_area_m2:.4f} m^2")
print(f"Daily evaporative loss: {burn_area_m2:.4f} m^2 * 2500 cc/m^2/day = {daily_evaporative_loss:.2f} cc/day")
print(f"Hourly evaporative loss: {daily_evaporative_loss:.2f} / 24 = {hourly_evaporative_loss:.2f} cc/hr\n")


# --- Step 3: Calculate Total Hourly Rate ---
total_hourly_rate = hourly_maintenance_fluid + hourly_evaporative_loss

print("Step 3: Calculate Total Maintenance Fluid Rate")
print("The total rate is the sum of the baseline rate and the evaporative loss rate.")
print(f"Total Rate (cc/hr) = {hourly_maintenance_fluid:.2f} + {hourly_evaporative_loss:.2f}")
print(f"Final Answer: The required maintenance fluid rate is {total_hourly_rate:.2f} cc/hr.")

final_answer = round(total_hourly_rate, 2)