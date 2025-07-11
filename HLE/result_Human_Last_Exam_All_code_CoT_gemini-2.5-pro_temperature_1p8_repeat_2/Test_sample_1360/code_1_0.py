import math

def calculate_time_resolution_requirement():
    """
    Calculates the average time between radioactive decay events for a given activity
    and explains why this is the dominant factor for time resolution.
    """
    # The given activity of the Bi 207 source is 1 kBq.
    activity_in_kBq = 1

    # Convert the activity from kilo-Becquerels (kBq) to Becquerels (Bq).
    # 1 kBq = 1000 Bq, and 1 Bq = 1 decay per second.
    activity_in_Bq = activity_in_kBq * 1000

    # The activity is the rate of decays (R). The average time interval (T)
    # between consecutive decays is the reciprocal of the rate.
    # T = 1 / R
    average_time_interval = 1 / activity_in_Bq

    # Convert the result to milliseconds for easier interpretation.
    # 1 second = 1000 milliseconds.
    average_time_interval_ms = average_time_interval * 1000

    print("Step 1: The activity of the source is given as {} kBq.".format(activity_in_kBq))
    print("Step 2: This activity corresponds to {} decays per second.".format(activity_in_Bq))
    print("Step 3: To resolve individual decay events, the measurement time resolution must be shorter than the average time between decays.")
    print("Step 4: We calculate the average time between decays as the reciprocal of the activity.")
    print(f"   Calculation: Average time = 1 / {activity_in_Bq} s^-1 = {average_time_interval:.4f} s")
    print("\nThis means, on average, a decay occurs every {:.1f} milliseconds.".format(average_time_interval_ms))
    print("\nConclusion: The rate of decay (the activity) dictates the time window one has to distinguish one event from the next.")
    print("This time (~1 ms) is much longer than other time scales in the experiment, like the electron's time-of-flight (~3-5 nanoseconds).")
    print("Therefore, the dominant factor that sets the minimum time resolution requirement is the measured activity of the source.")

calculate_time_resolution_requirement()