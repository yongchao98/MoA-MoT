def calculate_timescale_from_activity(activity_in_kBq):
    """
    Calculates the average time between decays for a given source activity
    and explains why activity is the dominant factor for time resolution.
    """
    
    # 1 kBq = 1000 Bq (decays per second)
    activity_in_Bq = activity_in_kBq * 1000
    
    # The average time between two independent decay events is the inverse of the activity.
    if activity_in_Bq > 0:
        average_time_s = 1 / activity_in_Bq
    else:
        average_time_s = float('inf')
        
    # Convert to milliseconds for readability
    average_time_ms = average_time_s * 1000
    
    # --- Output ---
    print(f"The source activity is {activity_in_kBq} kBq, which equals {activity_in_Bq} decays per second.")
    print("\nThe minimum time resolution requirement is dominantly set by the rate of decay events.")
    print("This rate can be used to find the average time between consecutive, independent decays.")
    print("\nCalculation:")
    print(f"  Average Time = 1 / Activity")
    print(f"  Average Time = 1 / {activity_in_Bq} decays/second")
    print(f"  Average Time = {average_time_s:.4f} seconds")
    print(f"  Average Time = {average_time_ms:.1f} milliseconds\n")
    
    print("Explanation:")
    print(f"To measure electrons from each decay 'individually', the detector system must be able to resolve events that are separated in time. "
          f"With an average time of {average_time_ms:.1f} ms between decays, the system's time resolution must be significantly shorter than this value to prevent unrelated events from being counted as one (pile-up).")
    print("Therefore, the measured activity of the source (1 kBq) is the dominant factor that determines this minimum requirement.")

# Run the calculation for the given problem
calculate_timescale_from_activity(1)