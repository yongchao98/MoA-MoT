import pint

def calculate_time_between_decays():
    """
    Calculates the average time between decays for a given activity.
    The dominant factor that sets the minimum time resolution requirement
    is the source activity, as it determines the average rate of events.
    To measure events individually, the system must be able to resolve them in time,
    meaning its resolution must be much smaller than the average time between events
    to keep the probability of "pile-up" low.
    """
    # Using pint for unit handling to ensure correctness
    ureg = pint.UnitRegistry()

    # Define the activity
    activity = 1 * ureg.kilobecquerel
    print(f"Source Activity: {activity:.0f}")

    # Convert to Becquerel (decays/second)
    activity_in_bq = activity.to(ureg.becquerel)
    print(f"This is equal to {activity_in_bq.magnitude:.0f} decays per second.")

    # Calculate the average time between decays
    # Average time = 1 / rate
    average_time = 1 / activity_in_bq
    print(f"\nThe average time between two consecutive decays is the reciprocal of the activity.")
    print(f"Average time = 1 / {activity_in_bq.magnitude:.0f} per second")

    # Convert to a more readable unit (milliseconds)
    average_time_ms = average_time.to(ureg.millisecond)
    print(f"Result: {average_time_ms.magnitude:.1f} milliseconds.")

    print("\nThis time (1 ms) is the characteristic timescale of the experiment.")
    print("The time resolution must be significantly better than this to individually measure electrons from different decays.")
    print("Therefore, the dominant factor setting this requirement is the activity of the source.")

calculate_time_between_decays()