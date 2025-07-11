import math

# --- Given Parameters ---
# D. The measured activity of the source (1 kBq)
activity_Bq = 1000.0  # 1 kBq in decays/second
# A. Distance between detectors
distance_detectors_m = 1.0  # 1 meter

# --- Other Physical Parameters ---
# Use a typical high-end energy for Bi-207 conversion electrons to calculate flight time
electron_kinetic_energy_MeV = 1.0
c_m_per_s = 299792458.0  # Speed of light in m/s
m_e_MeV_per_c2 = 0.511  # Electron rest mass in MeV/c^2

# Step 1: Calculate the timescale determined by the source activity.
# To resolve individual decays, the time resolution must be better than the
# average time separating them.
# Equation: T_avg_decay = 1 / Activity
average_time_between_decays_s = 1.0 / activity_Bq

print("--- Factor 1: Source Activity ---")
print(f"The source activity is {activity_Bq} decays/second.")
print("The average time between two consecutive decays is calculated as:")
print(f"1 / {int(activity_Bq)} decays/s = {average_time_between_decays_s} s")
print(f"This is equal to {average_time_between_decays_s * 1e3} ms (milliseconds).\n")

# Step 2: Calculate the timescale determined by the detector distance.
# This is the time it takes for a particle to travel from the source to a detector.
distance_source_to_detector_m = distance_detectors_m / 2.0

# Calculate electron speed (relativistically)
total_energy_MeV = electron_kinetic_energy_MeV + m_e_MeV_per_c2
gamma = total_energy_MeV / m_e_MeV_per_c2
beta_squared = 1.0 - (1.0 / (gamma**2))
beta = math.sqrt(beta_squared) # beta = v/c
electron_speed_m_per_s = beta * c_m_per_s
time_of_flight_s = distance_source_to_detector_m / electron_speed_m_per_s

print("--- Factor 2: Detector Distance ---")
print(f"The distance from the source to one detector is {distance_source_to_detector_m} m.")
print(f"The time of flight for a {electron_kinetic_energy_MeV} MeV electron is calculated as:")
print(f"Distance / Speed = {distance_source_to_detector_m} m / ({beta:.3f} * c) = {time_of_flight_s:.3e} s")
print(f"This is equal to {time_of_flight_s * 1e9:.3f} ns (nanoseconds).\n")


# Step 3: Compare and conclude.
print("--- Conclusion ---")
print(f"Timescale from Activity: {average_time_between_decays_s * 1e3:,.1f} ms")
print(f"Timescale from Distance: {time_of_flight_s * 1e9:.3f} ns")
print("\nThe time between decay events (milliseconds) is roughly 500,000 times longer than the electron flight time (nanoseconds).")
print("To individually measure electrons from different decays, the system must be able to distinguish between two events that occur randomly in time. The rate of these random events is defined by the source activity.")
print("A higher activity would lead to a shorter average time between events, thus demanding better time resolution to avoid 'pile-up'. The flight time is merely a small delay. Therefore, the activity is the dominant factor.")
<<<D>>>