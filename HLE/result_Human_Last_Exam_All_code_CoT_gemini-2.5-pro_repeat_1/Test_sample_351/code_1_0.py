import math

# Patient and burn details
weight_kg = 25
tbsa_percent = 45
bsa_m2 = 0.95

# Constants for the Galveston formula
MAINTENANCE_RATE_PER_M2_BSA = 1500  # mL/m^2/day
EVAPORATIVE_LOSS_RATE_PER_M2_BURN = 4000 # mL/m^2/day
HOURS_IN_DAY = 24

# --- Calculations ---

# 1. Daily maintenance fluid
maintenance_fluid_per_day = MAINTENANCE_RATE_PER_M2_BSA * bsa_m2

# 2. Daily evaporative loss fluid
tbsa_decimal = tbsa_percent / 100
burned_area_m2 = bsa_m2 * tbsa_decimal
evaporative_loss_per_day = EVAPORATIVE_LOSS_RATE_PER_M2_BURN * burned_area_m2

# 3. Total daily fluid
total_fluid_per_day = maintenance_fluid_per_day + evaporative_loss_per_day

# 4. Hourly rate
hourly_rate = total_fluid_per_day / HOURS_IN_DAY

# --- Output the results ---
print("Calculating Maintenance Fluid Needs for a Pediatric Burn Patient")
print("-------------------------------------------------------------")
print(f"Patient BSA: {bsa_m2} m^2")
print(f"TBSA Burn: {tbsa_percent}%")
print("\nThe calculation combines baseline maintenance needs with evaporative losses from the burn.")

print("\nStep 1: Calculate daily maintenance fluids (1500 mL/m^2 * BSA)")
print(f"Calculation: {MAINTENANCE_RATE_PER_M2_BSA} mL/m^2 * {bsa_m2} m^2 = {maintenance_fluid_per_day:.2f} mL/day")

print("\nStep 2: Calculate daily evaporative loss (4000 mL/m^2 * Burned Area)")
print(f"Burned Area = {bsa_m2} m^2 * {tbsa_decimal} = {burned_area_m2:.4f} m^2")
print(f"Calculation: {EVAPORATIVE_LOSS_RATE_PER_M2_BURN} mL/m^2 * {burned_area_m2:.4f} m^2 = {evaporative_loss_per_day:.2f} mL/day")

print("\nStep 3: Calculate total daily fluid requirement")
print(f"Calculation: {maintenance_fluid_per_day:.2f} mL/day + {evaporative_loss_per_day:.2f} mL/day = {total_fluid_per_day:.2f} mL/day")

print("\nStep 4: Convert to an hourly rate (cc/hr)")
print("Final Equation: ( (Maintenance Fluid) + (Evaporative Loss) ) / 24 hours")
# Using the numbers from the problem description in the final equation
print(f"Equation: ( ( {MAINTENANCE_RATE_PER_M2_BSA} * {bsa_m2} ) + ( {EVAPORATIVE_LOSS_RATE_PER_M2_BURN} * ( {bsa_m2} * {tbsa_percent/100.0} ) ) ) / {HOURS_IN_DAY}")
print(f"Result: {total_fluid_per_day:.2f} / {HOURS_IN_DAY} = {hourly_rate:.2f} cc/hr")

print(f"\nThe required maintenance fluid rate is approximately {math.ceil(hourly_rate)} cc/hr (rounded up).")