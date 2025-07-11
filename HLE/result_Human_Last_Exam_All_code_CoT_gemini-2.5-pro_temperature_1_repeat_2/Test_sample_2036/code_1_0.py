import math

# Step 1: Define assumptions based on standard ATC procedures and aircraft performance.
# These values are critical for the calculation.
SEPARATION_MINIMA_NM = 3.0  # Standard radar separation in nautical miles
APPROACH_SPEED_KTS = 140.0  # Assumed speed for the arriving aircraft in knots
DEPARTURE_AVG_SPEED_KTS = 160.0 # Assumed average speed for the departing aircraft during initial climb
TAKEOFF_PREP_TIME_SEC = 30.0    # Time from takeoff clearance to starting the roll
TAKEOFF_ROLL_TIME_SEC = 40.0    # Time for the takeoff roll until airborne

# Step 2: Calculate the total time needed for the departure to be safely separated.

# Time for the departing aircraft to fly the required 3 NM separation distance after getting airborne.
# Time (hours) = Distance / Speed
time_to_fly_3nm_hr = SEPARATION_MINIMA_NM / DEPARTURE_AVG_SPEED_KTS
# Convert this time to seconds for easier addition.
time_to_fly_3nm_sec = time_to_fly_3nm_hr * 3600

# Total time from giving takeoff clearance until the departing aircraft is 3 NM away.
total_time_for_departure_sec = TAKEOFF_PREP_TIME_SEC + TAKEOFF_ROLL_TIME_SEC + time_to_fly_3nm_sec

# Step 3: Calculate how far the arriving aircraft travels in that total time.
# This distance is the required separation from the VOR at the moment of takeoff clearance.

# Convert the total time back to hours to use with the speed in knots (NM/hour).
total_time_for_departure_hr = total_time_for_departure_sec / 3600

# Distance (NM) = Speed (knots) * Time (hours)
required_distance_from_vor_nm = APPROACH_SPEED_KTS * total_time_for_departure_hr

# Final Output
# The request is to show the final equation with the numbers.
# The calculation is: Required Distance = Arrival_Speed * (Prep_Time + Roll_Time + Flight_Time_for_Separation)
# Let's print the explanation and the result.

print("To ensure a 3 NM separation, we must calculate the total time for the departure and see how far the arrival travels in that time.")
print(f"Assumptions:")
print(f"- Arriving aircraft speed: {APPROACH_SPEED_KTS} kts")
print(f"- Departing aircraft average climb speed: {DEPARTURE_AVG_SPEED_KTS} kts")
print(f"- Required separation: {SEPARATION_MINIMA_NM} NM")
print(f"- Time from clearance to takeoff roll: {TAKEOFF_PREP_TIME_SEC}s")
print(f"- Time of takeoff roll: {TAKEOFF_ROLL_TIME_SEC}s\n")

print("Calculation:")
# We reconstruct the equation in the print statement for clarity, as requested.
# The core math is: required_distance = 140 * ((30 + 40 + (3 / 160 * 3600)) / 3600)
print(f"Required distance = {APPROACH_SPEED_KTS} kts * (({TAKEOFF_PREP_TIME_SEC}s + {TAKEOFF_ROLL_TIME_SEC}s + {round(time_to_fly_3nm_sec, 1)}s) / 3600s/hr)")
print(f"Required distance = {APPROACH_SPEED_KTS} kts * ({round(total_time_for_departure_sec, 1)}s / 3600s/hr)")
print(f"The arriving traffic must be at least {round(required_distance_from_vor_nm, 1)} NM from the VOR to clear the departing traffic.\n")
print(f"Final Answer: For separation, the arriving traffic needs to be {round(required_distance_from_vor_nm, 1)} miles from the VOR.")

# For the final answer format
final_answer = round(required_distance_from_vor_nm, 1)
# The prompt doesn't ask for this to be printed, but it's good practice to have it as a variable.
# print(f'<<<{final_answer}>>>') is not needed here as the format is handled outside the code block.