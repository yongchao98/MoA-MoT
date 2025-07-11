def calculate_time_resolution_factor(activity_kBq):
    """
    Calculates the average time between decays for a given activity.
    This demonstrates the dominant factor setting the time resolution requirement.
    """
    # Convert kilo-Becquerel (kBq) to Becquerel (Bq, decays per second)
    activity_Bq = activity_kBq * 1000

    # The rate of decays is the activity in Bq
    decay_rate = activity_Bq  # in decays/second

    # The average time between two consecutive decay events is the inverse of the rate
    average_time_s = 1 / decay_rate

    # Convert seconds to milliseconds for easier interpretation
    average_time_ms = average_time_s * 1000

    print("The activity of a radioactive source tells us the average number of decays per second.")
    print(f"An activity of {activity_kBq} kBq means {activity_Bq} * 1000 = {activity_Bq * 1000} decays per second on average.")
    print("\nThe minimum time resolution required to distinguish individual decay events is primarily determined by how frequently these events occur.")
    print("We can calculate the average time between consecutive decays:")
    print("\nEquation:")
    print(f"Average Time = 1 / Activity_in_Bq")
    print(f"\nResult:")
    # The final print statement shows the equation with the specific numbers
    print(f"Average Time = 1 / {activity_Bq} = {average_time_s:.3f} seconds, or {average_time_ms:.0f} milliseconds.")
    print("\nTo resolve most events individually, the detector's time resolution must be significantly smaller than this average time. Therefore, the source activity is the dominant factor.")


# Run the calculation for the given problem
calculate_time_resolution_factor(1)
