import math

def calculate_time_resolution_factor():
    """
    Calculates the average time between decays for a radioactive source
    to demonstrate why activity is the dominant factor for time resolution.
    """
    # The activity of the Bi 207 source is 1 kBq.
    activity_kBq = 1.0

    # Convert activity to Bq (decays per second).
    # 1 kBq = 1000 Bq
    activity_Bq = activity_kBq * 1000

    # The activity (A) is the number of decays per unit time.
    # The average time (T_avg) between two consecutive random decays is the
    # reciprocal of the activity.
    # This T_avg is the characteristic time that the measurement system
    # must be able to resolve to distinguish separate decay events.
    average_time_s = 1.0 / activity_Bq

    # Convert the time to more intuitive units like milliseconds (ms) and nanoseconds (ns)
    average_time_ms = average_time_s * 1000
    average_time_ns = average_time_s * 1e9

    print("Step 1: Define the source activity in decays per second (Bq).")
    print(f"Source Activity = {activity_kBq} kBq = {int(activity_Bq * 1000)} Bq\n")

    print("Step 2: Calculate the average time between consecutive decay events.")
    print("This is the dominant timescale for distinguishing separate events.")
    print("Formula: T_avg = 1 / Activity\n")

    print("Step 3: Perform the calculation.")
    print(f"T_avg = 1 / {int(activity_Bq)} s^-1")
    print(f"T_avg = {average_time_s} seconds")
    print(f"T_avg = {average_time_ms} milliseconds")
    print(f"T_avg = {int(average_time_ns)} nanoseconds\n")

    print("Conclusion:")
    print("The average time between decays is ~1 millisecond (1,000,000 ns).")
    print("Other factors like electron travel time are on the order of nanoseconds.")
    print("Since the events are separated by such a large average time, the rate of events (the activity) is the dominant factor setting the time resolution requirement.")

calculate_time_resolution_factor()