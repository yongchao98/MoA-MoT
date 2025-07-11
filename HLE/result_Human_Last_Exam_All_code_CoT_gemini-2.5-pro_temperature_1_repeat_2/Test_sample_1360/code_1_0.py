import math

def analyze_time_resolution_factors():
    """
    Analyzes the factors setting the time resolution requirement for detecting
    electrons from Bi-207 decay.
    """
    print("Analyzing the timescales involved in the experiment...")
    print("="*50)

    # 1. Timescale from Source Activity
    # The activity is the rate of decays. The average time between decays is its reciprocal.
    # This determines how frequently, on average, we must be able to record separate events.
    activity_kBq = 1.0
    activity_Bq = activity_kBq * 1000  # 1 kBq = 1000 decays/second

    print("Factor 1: Source Activity")
    print("The minimum time resolution must be sufficient to distinguish between separate decay events.")
    print("The average rate of these events is given by the source activity.")
    
    # Equation and Calculation
    avg_time_between_decays_s = 1 / activity_Bq
    print("\nCalculation of the average time between decays:")
    print(f"  Average Time = 1 / Activity")
    print(f"               = 1 / {activity_Bq} decays/s")
    print(f"               = {avg_time_between_decays_s} s")
    print(f"               = {avg_time_between_decays_s * 1000:.0f} milliseconds (ms)")

    print("-"*50)

    # 2. Timescale from Detector Distance
    # This determines the time-of-flight (ToF) of the electron from the source to the detector.
    # Bi-207 emits high-energy electrons that travel near the speed of light.
    distance_m = 1.0
    distance_to_detector_m = distance_m / 2
    speed_of_light_m_per_s = 3.0e8 # Approximate speed of high-energy electrons

    print("Factor 2: Distance between Detectors")
    print("The distance determines the electron's time-of-flight from the source to the detector.")
    
    # Equation and Calculation
    time_of_flight_s = distance_to_detector_m / speed_of_light_m_per_s
    print("\nCalculation of the time-of-flight:")
    print(f"  Time-of-Flight = Distance / Speed")
    print(f"                 = {distance_to_detector_m} m / {speed_of_light_m_per_s:.1e} m/s")
    print(f"                 = {time_of_flight_s:.2e} s")
    print(f"                 = {time_of_flight_s * 1e9:.2f} nanoseconds (ns)")
    
    print("="*50)

    # 3. Comparison and Conclusion
    print("Conclusion:")
    print(f"The average time between decay events is ~{avg_time_between_decays_s * 1000:.0f} ms.")
    print(f"The time it takes for an electron to reach a detector is ~{time_of_flight_s * 1e9:.2f} ns.")
    print("\nTo individually measure electrons from different decays, the system must be able to resolve events that are separated in time.")
    print("The random nature of radioactive decay means that while the *average* time between events is 1 ms, some events can occur much closer together.")
    print("The need to resolve these randomly occurring events sets the requirement for the time resolution. A higher activity would mean more events per second, a shorter average time between them, and thus a need for a better (smaller) time resolution to avoid event 'pile-up'.")
    print("Since the time between decays (1,000,000 ns) is about six orders of magnitude larger than the time-of-flight (~1.7 ns), the dominant factor setting the time resolution requirement is the rate of events, which is the measured activity of the source.")

analyze_time_resolution_factors()
<<<D>>>