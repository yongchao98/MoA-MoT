import math

# The activity is given in kilo-Becquerel (kBq).
activity_kBq = 1

# Convert activity to Becquerel (decays per second).
# 1 kBq = 1000 Bq
activity_Bq = activity_kBq * 1000

# The average time interval (T) between decays is the reciprocal of the activity (A).
# T = 1 / A
average_time_interval_s = 1 / activity_Bq

# For comparison, let's calculate the time of flight (ToF) for an electron.
# The distance from the source (in the middle) to a detector is 0.5 m.
distance_m = 0.5 
# The speed of an electron is a fraction of the speed of light. Let's use the speed of light as an upper bound for speed.
speed_of_light_m_per_s = 299792458
time_of_flight_s = distance_m / speed_of_light_m_per_s

print(f"The activity of the source is {activity_Bq} decays per second.")
print("The average time between two consecutive decay events is calculated as the reciprocal of the activity.")
print(f"Average time between decays = 1 / {activity_Bq} Bq = {average_time_interval_s} seconds")
print(f"This is equal to {average_time_interval_s * 1000} milliseconds.")
print("\nFor comparison:")
print(f"The time for an electron to travel 0.5 m to a detector is approximately {time_of_flight_s * 1e9:.2f} nanoseconds.")
print("\nConclusion:")
print("The average time between decays (~1 millisecond) is much longer than the electron's time of flight (~1.7 nanoseconds).")
print("Therefore, the rate of decay events (the activity) is the dominant factor that determines how fast the detection system needs to be to distinguish one event from the next.")
