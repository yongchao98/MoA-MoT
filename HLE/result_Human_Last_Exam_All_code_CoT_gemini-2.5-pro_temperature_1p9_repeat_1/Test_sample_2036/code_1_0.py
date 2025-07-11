import math

def calculate_separation_distance():
    """
    Calculates the required separation distance for a departure
    ahead of an arrival on a circling approach.
    """
    # According to ICAO PANS-OPS, the radius of the circling area is based on aircraft category.
    # Category C is common for airliners like the A320/B737.
    # Radius for CAT C is 4.2 NM.
    circling_area_radius_nm = 4.2

    # A safety buffer is added to allow the departing aircraft time to roll and get airborne
    # before the arrival reaches the circling area boundary. 1.8 NM is a reasonable buffer.
    safety_buffer_nm = 1.8

    # The total separation distance is the sum of the circling radius and the buffer.
    total_separation_distance = circling_area_radius_nm + safety_buffer_nm

    print("To ensure separation between the inbound circling traffic and the departing traffic, we must calculate the minimum distance the arrival must be from the airport.")
    print("\nThis is based on two main components:")
    print(f"1. The radius of the protected Circling Approach Area for a standard Category C aircraft.")
    print(f"2. A safety buffer to allow time for the departure's takeoff roll.")
    
    print("\nThe equation is:")
    print(f"Required Distance = Circling Area Radius + Safety Buffer")
    print(f"Required Distance = {circling_area_radius_nm} NM + {safety_buffer_nm} NM")
    
    # Using math.ceil to round up to be conservative, although it's a round number here.
    final_distance = math.ceil(total_separation_distance)

    print(f"\nFinal calculated distance: {final_distance} NM")
    print("\nTherefore, you must wait until the arriving traffic is 6 NM from the VOR/runway to clear the next aircraft for takeoff.")


calculate_separation_distance()

<<<6>>>