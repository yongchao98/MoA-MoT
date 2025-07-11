import math

# --- Constants ---
ACTIVITY_Bq = 1000  # 1 kBq = 1000 decays/second
DISTANCE_m = 0.5  # Distance from source to one detector
SPEED_OF_LIGHT_m_s = 3.0e8

# --- Electron Speed Calculation (for a 1 MeV electron as an example) ---
# Relativistic calculation: E_kinetic = (gamma - 1) * m_e * c^2
# We can estimate the speed to be a significant fraction of c without a full calculation,
# but for completeness, let's use a known value.
# A 1 MeV electron has a speed of approximately 0.94 times the speed of light.
ELECTRON_SPEED_fraction_of_c = 0.94
electron_speed_m_s = ELECTRON_SPEED_fraction_of_c * SPEED_OF_LIGHT_m_s

# --- Time Scale Calculations ---

# 1. Time scale determined by the source activity
# This is the average time between consecutive decay events.
avg_time_between_decays_s = 1.0 / ACTIVITY_Bq
avg_time_between_decays_ms = avg_time_between_decays_s * 1000

# 2. Time scale determined by the distance to the detectors
# This is the time it takes for an emitted electron to travel to a detector.
time_of_flight_s = DISTANCE_m / electron_speed_m_s
time_of_flight_ns = time_of_flight_s * 1e9

# --- Output Results and Comparison ---
print("Analysis of Time Scales:")
print("-" * 30)

print(f"Source Activity: {ACTIVITY_Bq} Bq")
print(f"Average time between decays = 1 / {ACTIVITY_Bq} Bq")
print(f"                            = {avg_time_between_decays_s:.4f} s")
print(f"                            = {avg_time_between_decays_ms:.1f} ms\n")


print(f"Distance to detector: {DISTANCE_m} m")
print(f"Electron Time-of-Flight = {DISTANCE_m} m / ({ELECTRON_SPEED_fraction_of_c:.2f} * c)")
print(f"                        = {time_of_flight_s:.2e} s")
print(f"                        = {time_of_flight_ns:.2f} ns\n")


# --- Conclusion ---
comparison_ratio = avg_time_between_decays_s / time_of_flight_s
print("Comparison:")
print(f"The average time between decays ({avg_time_between_decays_ms:.1f} ms) is ~{int(comparison_ratio):,} times larger than the electron time-of-flight ({time_of_flight_ns:.2f} ns).")
print("\nConclusion: To distinguish individual decay events from each other, the detector's time resolution must be able to handle events separated by an average of 1 ms. This is the dominant limiting factor set by the source's activity. The time-of-flight is a much shorter time scale and thus a less demanding constraint.")
