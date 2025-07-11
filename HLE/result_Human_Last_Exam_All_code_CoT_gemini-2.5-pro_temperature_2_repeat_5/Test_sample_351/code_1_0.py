import math

# Patient Information
# BSA (Body Surface Area) in square meters
bsa = 0.95

# Constants for Calculation
# Standard rate for maintenance fluids using the BSA method is 1500 mL/m^2/day
maintenance_rate_per_bsa_per_day = 1500
# Number of hours in a day
hours_in_day = 24

# --- Calculation ---
# 1. Calculate the total daily maintenance fluid volume
total_daily_fluid = bsa * maintenance_rate_per_bsa_per_day

# 2. Calculate the hourly rate by dividing the daily volume by 24
hourly_fluid_rate = total_daily_fluid / hours_in_day

# --- Output ---
print("This script calculates pediatric maintenance fluid needs using the Body Surface Area (BSA) method.")
print("\n--- Calculation Steps ---")
print(f"1. Daily Fluid Needs = BSA * Rate per Day")
print(f"   Daily Fluid Needs = {bsa} m^2 * {maintenance_rate_per_bsa_per_day} mL/m^2/day = {total_daily_fluid:.2f} mL/day")
print("\n2. Hourly Fluid Rate = Daily Fluid Needs / 24 hours")
# Using the variables in the final print statement as requested
print(f"   Hourly Fluid Rate = ({bsa} * {maintenance_rate_per_bsa_per_day}) / {hours_in_day}")

# Final Answer
# Using math.ceil to round up to the nearest tenth for clinical practicality.
final_answer = math.ceil(hourly_fluid_rate * 10) / 10
print(f"\nThe required maintenance fluid rate is: {final_answer} cc/hr")

# The final result in the requested format will be rounded to one decimal place.