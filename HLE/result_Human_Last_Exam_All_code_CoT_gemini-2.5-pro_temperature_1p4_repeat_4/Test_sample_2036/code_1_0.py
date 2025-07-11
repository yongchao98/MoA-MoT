import math

def calculate_separation_distance():
    """
    Calculates the required separation distance for an arriving aircraft
    to allow a departure from the opposite runway.
    """

    # --- Assumptions based on standard aviation practice ---
    # The typical speed of an airliner on final approach in knots (nautical miles per hour).
    inbound_speed_kts = 150
    
    # The time required for a departing aircraft from takeoff clearance to being safely airborne, in minutes.
    # This includes line-up and takeoff roll. 1.5 minutes (90 seconds) is a safe, conservative estimate.
    takeoff_time_min = 1.5
    
    # The standard safety buffer distance from the runway threshold that the arriving aircraft
    # must maintain until the departure is rolling. This is a common ATC "gate".
    final_separation_gate_nm = 4.0

    # --- Calculations ---
    
    # Convert takeoff time from minutes to seconds
    takeoff_time_sec = takeoff_time_min * 60
    
    # Calculate the distance the inbound aircraft travels during the departure's takeoff time.
    # Speed (NM/sec) = Speed (NM/hr) / 3600 (sec/hr)
    # Distance (NM) = Speed (NM/sec) * Time (sec)
    distance_traveled_nm = (inbound_speed_kts / 3600) * takeoff_time_sec
    
    # The total required distance is the distance the inbound travels plus the final safety gate.
    total_distance_nm = distance_traveled_nm + final_separation_gate_nm

    # --- Output Results ---
    print("--- Separation Calculation for Bilbao Approach ---")
    print(f"Assumed inbound aircraft speed: {inbound_speed_kts} kts")
    print(f"Assumed time for departure takeoff: {takeoff_time_min} minutes ({int(takeoff_time_sec)} seconds)")
    print(f"Required final separation gate: {final_separation_gate_nm} NM")
    print("\nStep 1: Calculate distance traveled by inbound during takeoff roll.")
    print(f"   Calculation: ({inbound_speed_kts} kts / 3600 s/hr) * {int(takeoff_time_sec)} s = {distance_traveled_nm:.2f} NM")
    
    print("\nStep 2: Add the safety buffer distance.")
    
    # The final print statement showing the full equation as requested.
    print(f"\nFinal Equation: (Distance Traveled) + (Safety Gate)")
    print(f"Result: {distance_traveled_nm:.2f} NM + {final_separation_gate_nm} NM = {total_distance_nm:.2f} NM")

    # Final rounded answer for clarity
    final_answer = round(total_distance_nm, 1)
    print(f"\nTherefore, the inbound traffic needs to be at least {final_answer} NM from the VOR.")
    print(f"<<<{final_answer}>>>")


if __name__ == "__main__":
    calculate_separation_distance()