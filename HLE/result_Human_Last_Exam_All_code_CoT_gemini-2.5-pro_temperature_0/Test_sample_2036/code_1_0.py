import math

# --- Step 1: Define assumptions and given data ---

# Time needed for the departing aircraft to take off and clear the runway path.
# This is a standard operational estimate for safety.
time_for_departure_to_clear_minutes = 2.0

# The arriving aircraft is on an instrument approach for RWY 30 (heading ~300 degrees).
# A typical airspeed for an aircraft on this initial approach segment.
arrival_airspeed_knots = 180.0

# The wind is from 120 degrees at 11 knots. This is a direct tailwind for an aircraft
# approaching RWY 30, so it increases the ground speed.
tailwind_knots = 11.0

# --- Step 2: Calculate the arriving aircraft's ground speed ---

# Ground Speed = Airspeed + Tailwind
arrival_ground_speed_knots = arrival_airspeed_knots + tailwind_knots

# --- Step 3: Calculate the required separation distance ---

# Convert the time from minutes to hours to match the speed unit (knots = nautical miles per hour).
time_for_departure_to_clear_hours = time_for_departure_to_clear_minutes / 60.0

# Use the formula: Distance = Speed * Time
required_separation_miles = arrival_ground_speed_knots * time_for_departure_to_clear_hours

# --- Step 4: Print the explanation and the final equation ---

print("To allow a departure from RWY 30, the arriving traffic needs to be far enough away to provide safe time separation.")
print(f"A safe time buffer for a takeoff is {time_for_departure_to_clear_minutes} minutes.")
print("\nFirst, we calculate the aircraft's ground speed:")
print(f"Ground Speed (knots) = Assumed Airspeed ({arrival_airspeed_knots} kt) + Tailwind ({tailwind_knots} kt) = {arrival_ground_speed_knots} kt")

print("\nNext, we calculate the distance required to provide the time buffer:")
print("Formula: Distance = Ground Speed * Time")
print("\nFinal Equation:")
# The f-string below prints the full equation with the numbers used.
print(f"Required Distance (NM) = {arrival_ground_speed_knots} knots * ({time_for_departure_to_clear_minutes} minutes / 60) hours = {round(required_separation_miles, 1)} NM")
print("\nTherefore, you need the arriving traffic to be at least 6.4 nautical miles from the VOR to clear the next traffic for takeoff.")

# The final answer is formatted and rounded to one decimal place.
final_answer = round(required_separation_miles, 1)
# print(f"<<<{final_answer}>>>") # This line is for the final output format, but the user wants the value.