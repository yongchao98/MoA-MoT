import math

def calculate_separation_distance():
    """
    Calculates the required distance from the VOR for an arriving aircraft
    to allow a departure, based on standard ATC separation rules.
    """
    # Step 1: Define key parameters based on standard operational values.
    # Minimum radar separation required between aircraft in a terminal area.
    min_radar_separation_nm = 5.0
    # A typical speed for a commercial jet on final approach.
    arrival_speed_kts = 180.0
    # A typical time from takeoff clearance until the aircraft is airborne and climbing.
    departure_takeoff_time_min = 1.0

    print("--- ATC Separation Calculation ---")
    print("This calculation determines the minimum distance an arriving aircraft must be from the airport")
    print("to safely allow another aircraft to depart.")
    print("\n Assumptions:")
    print(f"- Minimum Radar Separation: {min_radar_separation_nm} NM")
    print(f"- Arriving Aircraft Speed: {arrival_speed_kts} knots")
    print(f"- Time for Departure to Get Airborne: {departure_takeoff_time_min} minute")
    print("-" * 34)

    # Step 2: Calculate how far the arriving aircraft travels during the departure's takeoff.
    # Convert speed from knots (NM per hour) to NM per minute.
    arrival_speed_nm_per_min = arrival_speed_kts / 60.0
    distance_covered_by_arrival_nm = arrival_speed_nm_per_min * departure_takeoff_time_min

    # Step 3: Calculate the total required distance.
    # This is the minimum separation plus the distance the arrival covers.
    total_required_distance_nm = min_radar_separation_nm + distance_covered_by_arrival_nm

    # Step 4: Print the final equation and the result.
    print("\nCalculation:")
    print("To ensure separation is maintained, we calculate the distance the arriving aircraft")
    print(f"will cover in the {int(departure_takeoff_time_min)} minute it takes for the other plane to depart.")
    print(f"This distance is then added to the required minimum separation of {int(min_radar_separation_nm)} NM.")

    print("\nFinal Equation:")
    print(f"Required Distance = Minimum Separation + (Arrival Speed / 60) * Departure Time")
    # Using integer casting for cleaner print output of the equation components
    print(f"Required Distance = {int(min_radar_separation_nm)} NM + ({int(arrival_speed_kts)} kts / 60) * {int(departure_takeoff_time_min)} min")
    print(f"Required Distance = {int(min_radar_separation_nm)} NM + {int(distance_covered_by_arrival_nm)} NM")
    print(f"Required Distance = {total_required_distance_nm} NM")

    print("\nConclusion:")
    print(f"The arriving traffic must be at least {total_required_distance_nm} NM from the VOR to clear the next traffic for takeoff.")
    
    return total_required_distance_nm

# Execute the calculation and store the final answer
final_answer = calculate_separation_distance()
print(f"\n<<<{final_answer}>>>")
